
      subroutine minime (ids,phct)
c-----------------------------------------------------------------------

      INCLUDE 'link_fnl_shared.h'

      USE LCONF_INT

      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer ids, i, j, k, l, ncon, ocon, nact, iact(2*m4),
     *        maxfcn, neq, nvar, ntot, phct, lda

      double precision lambda(10), plb(10), gsol1, gval, sum,
     *                 pub(10), gfinal, dbz, acc, ppp(10)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)


      integer jds, tphct
      double precision az, bz
      common/ cstexp /az(36,10),bz(36),jds,tphct

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct


      external gsol2

      external gsol1 
c-----------------------------------------------------------------------
c                                 create constraint matrix z
c                                 --------------------------
      ntot = nstot(ids)
      nvar = ntot - 1
c                                 initialize bounds
      do i = 1, nvar
         pub(i) = 1d99
         plb(i) = -1d99
      end do 

      if (.not.mus) then 
         write (*,*) 'no mus'
         call errpau
      end if 

      jds = ids
      tphct = phct
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

      write (*,*) sum, g2(phct)
      
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
c                                 tolerance
      acc = 1d-5
c                                 no linear constraint
      neq = 0
c                                 ncon - neq linear inequalities,
c                                 to be counted:
      ncon = neq
c                                 for each site
      do i = 1, msite(ids)
c                                 for each species
         do j = 1, zsp(ids,i) + 1
c                                 initial az, bz
            ncon = ncon + 1
            dbz = 0d0
c                                 both Temkin and non-Temkin have
c                                 -Az*p <= 0 constraints:
            bz(ncon) = dcoef(0,j,i,ids) - acc

            do k = 1, nvar

               az(ncon,k) = 0d0

            end do

            do k = 1, lterm(j,i,ids)

               if (ksub(k,j,i,ids).ne.ntot) then 

                  az(ncon,ksub(k,j,i,ids)) = az(ncon,k) 
     *                                     - dcoef(k,j,i,ids)

               else 
c                                 decompose p(ntot) into 1 - p(1) - ...- p(nvar)
                  dbz = dcoef(k,j,i,ids)

                  do l = 1, nvar

                     az(ncon,l) = az(ncon,l) + dbz

                  end do

               end if 

            end do

            bz(ncon) = bz(ncon) + dbz
c                                 non-Temkin have the Az*p <= 1 constraint
            if (zmult(ids,i).ne.0d0) then

               ocon = ncon

               ncon = ncon + 1

               bz(ncon) =  1d0 - dcoef(0,j,i,ids) - dbz - acc

               do k = 1, nvar

                  az(ncon,k) = -az(ocon,k)

               end do

            end if

         end do 

      end do

      maxfcn = 400

c lconf f77
c     CALL LCONF (gsol2, NVAR, NCON, NEQ, Az, LDA, Bz, pLB, pUB,
c    *         ppp, ACC, MAXFCN, ppp, gfinal, NACT, IACT, lambda)

c lconf f90:
      call lconf (gsol2,neq,az,bz,plb,pub,ppp,xguess=ppp,
     *            maxfcn = maxfcn, acc = acc, obj = gfinal,
     *            nact = nact, iact = iact, alamda = lambda)

      end 

      subroutine gsol2 (nstt,ppp,gval)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, nstt

      logical bad

      double precision ppp(*), gval, gsol1, psum

      external gsol1 

      integer jds, tphct
      double precision az, bz
      common/ cstexp /az(36,10),bz(36),jds,tphct

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

      do i = 1, nstot(jds) - 1
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

1000  format (g12.6,1x,12(f7.5,1x))
1010  format (2(g14.7,2x))

      end



