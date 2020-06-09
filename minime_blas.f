      subroutine dummy
      end 

      subroutine minfrc (ids,kds)
c-----------------------------------------------------------------------
c minimize the omega function for the independent endmember fractions
c of solution ids subject to site fraction constraints

c     number of independent endmember fractions -> nstot-1 (<m19)
c     number of independent endmember fractions -> nz (<m20)
c     closure is forced in the objective function (gsol2)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, kds, nvar, iter, iwork(m22),
     *        istuff(10), istate(m21), istart, nclin, ntot
c DEBUG691
     *         , i, j, iprint

      double precision ggrd(m19), lapz(m20,m19),
     *                 bl(m21), bu(m21), gfinal, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),stuff(1)
c DEBUG691                    dummies for NCNLN > 0
     *                 ,c(1),cjac(1,1),yt(m14),zp,zt(m10,m11)

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      external gsol2, dummy

      data iprint/0/

      save iprint
c-----------------------------------------------------------------------
      if (.not.mus) then 
         write (*,*) 'no mus'
         call errpau
      end if

      nclin = nz(ids)
      ntot = nstot(ids)
      nvar = ntot - 1
c                                 finite difference increments
c                                 will be estimated at this 
c                                 coordinate, so choose a feasible 
c                                 composition
      ppp(1:nvar) = pa(1:nvar)
c                                 initialize bounds
      bu(1:nvar) = 1d20
      bl(1:nvar) = -1d20
c                                 load the local constraints 
c                                 from the global arrays
      lapz(1:nclin,1:nvar) = apz(ids,1:nclin,1:nvar)

      bl(nvar+1:nvar+nclin) = zl(ids,1:nclin)
      bu(nvar+1:nvar+nclin) = zu(ids,1:nclin)

      istart = -1
c                                 the solution model index
      istuff(1) = ids
c                                 print/save obj value 
c     istuff(2) set by nlpopt
c                                 obj call counter
      istuff(3) = 0
c                                 saved obj value counter
      istuff(4) = 0
c                                 refinement point index
      istuff(5) = kds

      CALL E04UEF ('nolist')
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

      call e04ucf (nvar,nclin,m20,m19,lapz,bl,bu,gsol2,
     *             iter,istate,clamda,gfinal,ggrd,r,ppp,iwork,
     *             m22,work,m23,istuff,stuff,istart,iprint)

c     call e04ucf (nvar,nclin,0,m20,1,m19,lapz,bl,bu,dummy,gsol2,
c    *             iter,istate,c,cjac,clamda,gfinal,ggrd,r,ppp,iwork,
c    *             m22,work,m23,istuff,stuff,istart)

c     write (*,1000) 'pas',pa(1:3)
c     write (*,1000) 'p0as',p0a(1:3)
1000  format (a10,10(g14.7,1x))

c     do i = 1, nz(ids)
c        zp = 0d0 
c        do j = 1, nvar
c           zp = zp + lapz(i,j) * ppp(j)
c        end do

c        write (*,1010) i, bl(i),zp,bu(i)

c     end do

c     call p2z (pa,zt,ids)

1010  format (i5,1x,10(g14.7,1x))
      end

      subroutine gsol2 (mode,nvar,ppp,gval,ggrd,istart,istuff,stuff)
c-----------------------------------------------------------------------
c function to evaluate gibbs energy of a solution for nlpopt. can call 
c either gsol1 that does o/d or gsol4 that does not, gsol1 seems to give
c better results presumably because it's using analytical gradients.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, jds, nvar, mode, istuff(*), istart

      double precision ppp(*), gval, gsol4, ggrd(*), stuff(*),
     *                 gsol1, g, sum, scp(k5), zt(m10,m11), sump

      external gsol4, gsol1

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      jds = istuff(1)

      sum = 0d0

      do i = 1, nvar
         sum = sum + ppp(i)
         p0a(i) = ppp(i)
      end do

      p0a(nstot(jds)) = 1d0 - sum
      pp(1:nstot(jds)) = p0a(1:nstot(jds))
      pa(1:nstot(jds)) = p0a(1:nstot(jds))

      call makepp (jds)

c     write (*,1000) 0, (pa(i),i=1,nstot(jds))

      g = gsol1 (jds)

c     write (*,1000) 1, (pa(i),i=1,nstot(jds))
c                                 get the bulk composition from pp
      call getscp (scp,sum,jds,jds,.false.)

c     write (*,1000) 0, (scp(i),i=1,icp)

      gval = g

      do i = 1, icp
         gval = gval - scp(i)*mu(i)
c         if (scp(i).lt.-1d2*zero) then
c            write (*,*) 'gsol2, large -composition',jds,i,scp(i)
c            call p2z (pa,zt,jds)
c            scp(i) = 0d0
c         end if
      end do
c                                  normalize for appearances
      gval = gval/sum

c     write (*,1000) gval, g

      istuff(3) = istuff(3) + 1

      if (istuff(2).ne.0) then
c                                 save the composition
         istuff(4) = istuff(4) + 1
c        write (*,1000) gval, g/sum, (pa(i),i=1,nstot(jds))
c                                 there is a small possibility 
c                                 that the composition of a phase
c                                 may move entirely into the 
c                                 constrained composition space
         if (sum.gt.zero) then
c                                 increment the counter
            jphct = jphct + 1
c                                 the solution model pointer
            jkp(jphct) = jds
c                                 the refinement point pointer
            hkp(jphct) = istuff(5)
c                                 save the normalized g
            g2(jphct) = g/sum
c                                 save the normalized bulk
            cp2(1:icomp,jphct) = scp(1:icomp)/sum
c                                 sum scp(1:icp)
            c2tot(jphct) = sum
c                                 save the endmember fractions
            icoz(jphct) = zcoct

            zco(zcoct+1:zcoct+nstot(jds)) = pa(1:nstot(jds))

            zcoct = zcoct + nstot(jds)

         end if

      end if

1000  format (2(g12.6,1x),12(f8.5,1x))
1010  format (2(g14.7,2x))

      end

      subroutine gsol3 (mode,nvar,ppp,gval,ggrd,istart,istuff,stuff)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jds, nvar, istart, mode, istuff(*)

      double precision ppp(*), gval, psum, ggrd(*), stuff(*)

      double precision gsol4, omega

      external gsol4, omega

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c-----------------------------------------------------------------------
      jds = istuff(1)

      psum = 0d0

      do i = 1, nvar
         psum = psum + ppp(i)
         p0a(i) = ppp(i)
      end do

      p0a(nstot(jds)) = 1d0 - psum

      do i = 1, nstot(jds)
         if (dabs(p0a(i)).lt.zero) p0a(i) = 0d0
         pp(i)  = p0a(i)
         pa(i)  = p0a(i)
      end do

      call makepp (jds)

      if (istuff(3).eq.0d0) then 
         gval = gsol4 (jds)
      else
         gval = - omega(jds,p0a)
      end if

c     write (*,1000) 0, (pa(i),i=1,nstot(jds))

c     write (*,1000) 0, (cp2(i,tphct),i=1,icp)
c     write (*,1000) gval, (pa(i),i=1,nstot(jds))

1000  format (g12.6,1x,12(f8.5,1x))
1010  format (2(g14.7,2x))

      end




      double precision function gsol4 (id)
c-----------------------------------------------------------------------
c gsol4 computes the total (excess+ideal) free energy of solution
c for a solution identified by index ids for the speciation input
c via cxt7, i.e., in contrast to gsol1 it does not compute the 
c speciation of o/d models.

c gsol4 assumes the endmember g's have been calculated by gall.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      double precision gg

      double precision omega, gfluid, gzero,
     *                 gex, gfesi, gfesic, gfecr1, gerk, ghybrid, gfes

      external omega, gfluid, gzero, gex, gfesi,
     *         gfesic, gfecr1, gerk, ghybrid, gfes

      integer jend
      common/ cxt23 /jend(h9,m4)

      double precision g
      common/ cst2 /g(k1)

      double precision enth
      common/ cxt35 /enth(j3)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c----------------------------------------------------------------------
      gg = 0d0

      if (specil(id)) then
c                                 special is reserved for special models 
c                                 that set standard flags (lorder and/or lrecip)
c                                 currently only Nastia's version of BCC/FCC Fe-Si-C Lacaze and Sundman
            gg =  gfesic (y(1),y(3),y(4),
     *                    g(jend(id,3)),g(jend(id,4)),
     *                    g(jend(id,5)),g(jend(id,6)),ksmod(id))

      else if (simple(id)) then
c                                 -------------------------------------
c                                 macroscopic formulation for normal solutions.
         call gdqf (id,gg,p0a)

         gg = gg - t * omega(id,p0a) + gex(id,p0a)
c                                 get mechanical mixture contribution
         do k = 1, lstot(id)
            gg = gg + p0a(k) * g(jend(id,2+k))
         end do

      else if (lorder(id)) then
c                                 dqf, excess, and entropic effects.
         call gdqf (id,gg,pp)

         gg = gg - t * omega(id,p0a) + gex(id,p0a)
c                                 enthalpic effect of forming the ordered
c                                 species:
         do k = 1, nord(id)
            gg = gg + p0a(lstot(id)+k)*enth(k)
         end do
c                                 plus the mechanical mixture of the 
c                                 disordered endmembers
         do k = 1, lstot(id)
            gg = gg + g(jend(id,2+k)) * pp(k)
         end do

      else if (ksmod(id).eq.0) then
c                                 ------------------------------------
c                                 internal fluid eos
         gg = gfluid(y(1))

         do k = 1, 2
            gg = gg + gzero(jnd(k))*y(k)
         end do

      else if (ksmod(id).eq.20) then
c                                 electrolytic solution, need to check
c                                 that this thing is getting the right
c                                 partial molar volumes.
         call slvnt1 (gg)

         call slvnt2 (gg)

      else if (ksmod(id).eq.26) then
c                                 ------------------------------------
c                                 andreas salt model
         call hcneos (gg,y(1),y(2),y(3))

         do k = 1, 3
            gg = gg + y(k) * g(jend(id,2+k))
         end do

      else if (ksmod(id).eq.29) then
c                                 -------------------------------------
c                                 BCC Fe-Si Lacaze and Sundman
         gg =  gfesi(y(1),g(jend(id,3)),g(jend(id,4)))

      else if (ksmod(id).eq.32) then
c                                 -------------------------------------
c                                 BCC Fe-Cr Andersson and Sundman
         gg =  gfecr1(y(1),g(jend(id,3)),g(jend(id,4)))

      else if (ksmod(id).eq.39) then
c                                 -------------------------------------
c                                 generic hybrid EoS
c                                 initialize pointer array
         do k = 1, nstot(id)
c                                 sum pure species g's
            gg = gg + g(jnd(k)) * y(k)

         end do
c                                 compute and add in activities
         gg = gg + ghybrid (y)

      else if (ksmod(id).eq.41) then
c                                 hybrid MRK ternary COH fluid
         call rkcoh6 (y(2),y(1),gg)

         do k = 1, 3
            gg = gg + g(jnd(k)) * y(k)
         end do

      else if (ksmod(id).eq.40) then
c                                 MRK silicate vapor
         do k = 1, nstot(id)
            gg = gg + gzero (jnd(k)) * y(k)
         end do

         gg = gg + gerk(y)

      else if (ksmod(id).eq.42) then
c                                 ------------------------------------
c                                 Fe-S fluid (Saxena & Eriksson 2015)
         gg =  gfes(y(2),g(jend(id,3)),g(jend(id,4)))

      else

         write (*,*) 'what the **** am i doing here?'
         call errpau

      end if

      gsol4 = gg

      end




      subroutine makayx (id)
c----------------------------------------------------------------------
c subroutine to construct the ayx matrices for the y to x conversion
c of each subcomposition.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ii, i, j, k, l, m, id

      double precision xt

      integer mx
      double precision ayx
      common/ csty2x /ayx(h9,h4,mst*msp,m4),mx(h9,h4)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      do ii = 1, pop1(id)
c                                 for subcomposition ii 
c                                 the Ax matrix in Ax*y = x
c                                 m = sum( ispg(1:istg) )
c                                 n = pvert(2) - pvert(1) + 1

         mx(id,ii) = 0

         do i = 1, istg(id,ii)
            mx(id,ii) = mx(id,ii) + ispg(id,ii,i)
         end do 

         do i = 1, ncoor(id)
            do j = 1, pvert(id,ii,2) - pvert(id,ii,2) + 1
               ayx(id,ii,i,j) = 0d0 
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

                     ayx(id,ii,i+m,j) = 1d0

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

         l = 1
         m = 1

         do i = 1, mx(id,ii)

            xt = 0d0
c                                 the algebra
            do k = pvert(id,ii,1), pvert(id,ii,2)

               j = k - pvert(id,ii,1) + 1

               xt = xt + ayx(id,ii,i,j)*y(k)

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

      subroutine makapz (id)
c----------------------------------------------------------------------
c subroutine to construct the apz matrix for the p' to independent z 
c limits, where p' is the first nstot-1 elements of p.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id, i, j, k, l, m, nvar, zct

      double precision dbz

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c---------------------------------------------------------
      nvar = nstot(id) - 1
c                                 to be counted:
      nz(id) = 0
c                                 for each site
      do i = 1, msite(id)
c                                 for each species
         if (zmult(id,i).eq.0d0.or.ksmod(id).eq.688) then
            zct = zsp1(id,i)
         else 
            zct = zsp(id,i)
         end if

         do j = 1, zct
c                                 initial az, bz
            nz(id) = nz(id) + 1

            dbz = 0d0
c                                 both Temkin and non-Temkin have
c                                 Az*p >= 0 constraints:
            apz(id,nz(id),1:nvar) = 0d0

            do k = 1, lterm(j,i,id)

               m = ksub(k,j,i,id)

               if (m.ne.nstot(id)) then 

                  apz(id,nz(id),m) = apz(id,nz(id),m) 
     *                               + dcoef(k,j,i,id)

               else 
c                                 decompose p(ntot) into 1 - p(1) - ...- p(nvar)
                  dbz = dcoef(k,j,i,id)

                  do l = 1, nvar

                     apz(id,nz(id),l) = apz(id,nz(id),l) - dbz

                  end do

               end if 

            end do

            zl(id,nz(id)) = - dbz
c                                 non-Temkin have the Az*p <= 1 constraint
            if (zmult(id,i).ne.0d0) then

               zu(id,nz(id)) =  1d0 - dbz

            else 

               zu(id,nz(id)) = 1d20

            end if

         end do

         if (zmult(id,i).gt.0d0.and.ksmod(id).ne.688) then 
c                                 a pre-688 model, need to make
c                                 a limit expression for the missing
c                                 site fraction

c                                 pointer to the previous constraints for
c                                 the site
            k = nz(id) - zsp(id,i) + 1

            nz(id) = nz(id) + 1
            dbz = 1d0
            apz(id,nz(id),1:nvar) = 0d0

            do j = 1, zsp(id,i)
               dbz = dbz + zl(id,k)
               do l = 1, nvar
                  apz(id,nz(id),l) = apz(id,nz(id),l) - apz(id,k,l)
               end do
               k = k + 1
            end do

            zl(id,nz(id)) =  -dbz
            zu(id,nz(id)) =  1d0 - dbz

         end if

      end do

      end


      subroutine sety2x (id,y)
c----------------------------------------------------------------------
c subroutine to convert independent disordered y to subcomposition
c x's, assumes y's are normalized.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ii, i, j, k, l, m, id

      double precision xt,y(m4)

      integer mx
      double precision ayx
      common/ csty2x /ayx(h9,h4,mst*msp,m4),mx(h9,h4)

      double precision z, pa, p0a, x, w, yt, wl, pp
      common/ cxt7 /yt(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
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

         l = 1
         m = 1

         do i = 1, mx(id,ii)

            xt = 0d0
c                                 the algebra
            do k = pvert(id,ii,1), pvert(id,ii,2)

               j = k - pvert(id,ii,1) + 1

               xt = xt + ayx(id,ii,i,j)*y(k)

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

      subroutine p2yx (id)
c-----------------------------------------------------------------------
c converts the independent endmember fractions to 0-1 bounded barycentric 
c coordinates:

c     number of bounded coordinates -> mstot (<m4)
c     number of independent fractions -> nstot (<m14)
c     number of linear constraints -> the number of independent
c        site fractions + closure (<m20)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw,lw,mvar,mcon

      parameter (mvar=m4, mcon=m20, liw=2*mvar+3, 
     *           lw=2*(mcon+1)**2+7*mvar+5*mcon)

      integer ncon, id, is(mvar+mcon), iw(liw), idead, istart

      double precision ax(mcon), clamda(mvar+mcon), wrk(lw), c(mvar),
     *                 a(mcon,mvar), b(mcon)

      double precision ayz
      common/ csty2z /ayz(h9,m20,m4)

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ldq,ldt,ncolt
      common/ be04nb /ldt,ncolt,ldq
c-----------------------------------------------------------------------
c                                 get the disordered p's
      call minfxc (id,.true.)
c                                 get the site fraction constraints
      call p2zind (pa,b,id)
c                                 dummy objective function coefficients
c                                 (only 1 feasible point)
      c(1:mstot(id)) = 1d0
c                                 load the ayz constraint matrix
      a(1:nz(id),1:mstot(id)) = ayz(id,1:nz(id),1:mstot(id))
c                                 add the closure constraint
      ncon = nz(id) + 1
      a(ncon,1:mstot(id)) = 1d0
      b(ncon) = 1d0
c                                 cold start
      istart = 0
      ldt = ncon
      ldq = ldt

      if (lopt(28)) call begtim (2)
c                                 optimize by nag
      call lpnag (mstot(id),ncon,a,mcon,b,c,is,y,ax,
     *            clamda,iw,liw,wrk,lw,idead,l6,istart)

      if (lopt(28)) call endtim (2,.true.,'p2y inversion')

      if (idead.gt.0) then
c                                 look for severe errors
         call lpwarn (idead,'LPOPT ')
         call errpau

      end if
c                                 reset ldt, ldq, istart for phase eq
      ldq = icp + 1
      ldt = ldq
      istart = 0
c                                 convert the y's to x's
      call sety2x (id,y)

      end

      subroutine p2zind (p,z,ids)
c----------------------------------------------------------------------
c subroutine to compute independent site fractions (or molar amounts 
c for temkin sites) and load them sequentially into the 1d array z
c
c non-temkin models:
c     zsp - number of independent site fractions
c     z(l) - molar site fraction of species j on site i
c temkin models:
c     zsp - number of species
c     z(l) - molar amount of species j on site i
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision p(*), z(*)

      integer i,j,k,l,ids

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
c                                 for each site
      l = 0

      do i = 1, msite(ids)
c                                 get site fractions
         do j = 1, zsp(ids,i)

            l = l + 1

            z(l) = dcoef(0,j,i,ids)
c                                 for each term:
            do k = 1, lterm(j,i,ids)

               z(l) = z(l) + dcoef(k,j,i,ids) * p(ksub(k,j,i,ids))

            end do

         end do

      end do

      end

      subroutine p2zall (y,z,ldz,ids)
c----------------------------------------------------------------------
c subroutine to compute all site fractions (or molar amounts for temkin
c sites) computed from equations
c
c non-temkin models:
c     zsp1 - number of site fractions
c     zsp - number of independent site fractions zsp-1
c     z(i,j) - molar site fraction of species j on site i
c temkin models:
c     zsp1 - number of species
c     zsp - zsp1
c     z(i,j) - molar amount of species j on site i
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision y(*), zt, z(ldz,*)

      integer i,j,k,ldz,ids

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
c                                 for each site
      do i = 1, msite(ids)

         zt = 0d0

         if (zmult(ids,i).ne.0d0.and.ksmod(ids).ne.688) then
c                                 get site fractions
            do j = 1, zsp(ids,i)

               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) +
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do

               zt = zt + z(i,j)

            end do

            z(i,j) = 1d0 - zt

         else if (zsp1(ids,i).gt.1) then
c                                 temkin or 688 model format, all species fractions are available
            do j = 1, zsp1(ids,i)
c                                 molar site population
               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) + 
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do

            end do

         end if

      end do

      end

      subroutine p2z (y,z,ids)
c----------------------------------------------------------------------
c subroutine to compute site fractions computed from equations entered by
c user for configurational entropy (macroscopic form). with range checks.
c
c non-temkin models:
c     z(i,j) - molar site fraction of species j on site i.
c temkin models:
c     z(i,j) - molar amount of species j on site i.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical badz, bad

      external badz

      double precision y(m4), zt, z(m10,m11)

      integer i,j,k,ids

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c----------------------------------------------------------------------
      bad = .false.
c                                 for each site
      do i = 1, msite(ids)

         zt = 0d0

         if (zmult(ids,i).ne.0d0.and.ksmod(ids).ne.688) then
c                                 get site fractions
            do j = 1, zsp(ids,i)

               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) +
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do

               bad = badz(z(i,j))

               if (bad) exit 

               zt = zt + z(i,j)

            end do

            if (bad) exit

            z(i,j) = 1d0 - zt

            bad = badz(z(i,j))

         else if (zsp1(ids,i).gt.1) then
c                                 temkin or 688 model format, all species fractions are available
            do j = 1, zsp1(ids,i)
c                                 molar site population
               z(i,j) = dcoef(0,j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)

                  z(i,j) = z(i,j) + 
     *                     dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))

               end do

               write (*,1001) i,j,z(i,j)
c                                 non-temkin (688)
               if (zmult(ids,i).gt.0d0.and.badz(z(i,j))) then 
c DEBUG691
                     call warn (72,
     *                       zt,i,'the expression for z('//
     *                       znames(ids,i,j)//') on '//znames(ids,i,0)//
     *                       ' in '//' is incorrect.')


               end if

               zt = zt + z(i,j)

            end do

c           write (*,1001) i,j,1d0-zt

            if (ksmod(ids).eq.688.and.zmult(ids,i).gt.0d0) then 
c                                 non-temkin, fractions must sum to 1
               if (dabs(zt-1d0).gt.1d4*zero) then

                  write (*,'(/,a,g14.6)') 'site fraction sum = ',zt

                  call warn (72,zt,i,
     *                       'site fractions on '//znames(ids,i,0)// 
     *                       ' in  do not sum to 1.')

               end if

            else if (zt.gt.0d0) then
c                                 temkin, if site exists, check fractions
               do j = 1, zsp(ids,i)

                  bad = badz(z(i,j)/zt)

                  if (bad) exit

               end do

            else if (zt.lt.-zero) then
c                                 negative site?
               bad = .true.

            end if

         end if

         if (bad) exit 

      end do

1000  format (/,'**error ver071** during testing of dependent endmember'
     *       ,' ',a,' the following invalid site fraction (z = ',g12.6,
     *        ')',/,'was found. The cause of this error may be either ',
     *       'the dependent endmember definition or invalid site',/,
     *       'fraction expressions for one or more of the independent ',
     *       'endmembers of ',a,/)

1001  format (i2,1x,i2,2x,g14.6)

      end


      subroutine makayz (id)
c----------------------------------------------------------------------
c subroutine to make the ayc matrix for ayz*y = z, z is the independent
c subset of the site fractions. must be called after makapz (for nz).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, id

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision ayz
      common/ csty2z /ayz(h9,m20,m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      integer jend
      common/ cxt23 /jend(h9,m4)
c----------------------------------------------------------------------
      ayz(id,:,:) = 0d0
c                                 independent endmembers:
      do k = 1, lstot(id)
c                                 endmember is in column knsp(k,id)
         pa(:) = 0d0
         pa(k) = 1d0
         call p2zind (pa,z,id)
         ayz(id,1:nz(id),knsp(k,id)) = z(1:nz(id))

      end do
c                                 dependent endmembers:
      do k = 1, ndep(id)
c                                 the index of the dependent endmember
c                                 in y is
         j = knsp(lstot(id)+k,id)

         do l = 1, ndph(k)
c                                 the dependent disordered endmember decomposes
c                                 to independent disordered endmember idep(k,l):
c                                 this is insanely inefficient, but who cares?
            pa(:) = 0d0
            pa(iy2p(idep(k,l))) = 1d0
            call p2zind (pa,z,id)

            do i = 1, nz(id)
               ayz(id,i,j) = ayz(id,i,j) + nu(k,l)*z(i)
            end do

         end do

      end do

      end

      subroutine makayc (id)
c----------------------------------------------------------------------
c subroutine to make the ayc matrix for ayc*y = c
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, m, n, id

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision ayc,cyc
      common/ csty2c /ayc(h9,k5,m4),cyc(h9,m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer jend
      common/ cxt23 /jend(h9,m4)
c----------------------------------------------------------------------
      ayc(id,:,:) = 0d0
      cyc(id,:) = 1d0
c                                 independent endmembers:
      do k = 1, lstot(id)
c                                 endmember is in column j
         j = knsp(k,id)

         do i = 1, icomp
            ayc(id,i,j) = cp(i,jend(id,2+k))
         end do

      end do
c                                 dependent endmembers:
      do k = 1, ndep(id)
c                                 the index of the dependent endmember
c                                 in y is
         j = knsp(lstot(id)+k,id)

         do l = 1, ndph(k)

            if (idep(k,l).le.mstot(id)) then
c                                 the dependent disordered endmember decomposes
c                                 to independent disordered endmember idep(k,l):
               do i = 1, icomp
                  ayc(id,i,j) = ayc(id,i,j) + 
     *                          nu(k,l)*cp(i,jend(id,2+iy2p(idep(k,l))))
               end do

            else
c                                 the dependent disordered endmember decomposes
c                                 to ordered endmember idep(j,l) - mstot, 
c                                 decompose the ordered endmember into chemically
c                                 equivalent disordered endmembers:
               m = idep(k,l) - mstot(id)
             
               do n = 1, nr(m)
c                                  the endmember decomposes to iddeps(n,m) in the
c                                  y array, or ideps(n,m,id) in the p array.
                  cyc(id,ideps(n,m,id)) = 0d0

                  do i = 1, icomp

                     ayc(id,i,j) = ayc(id,i,j) + nu(k,l)*depnu(n,m)
     *                           * cp(i,jend(id,2+ideps(n,m,id)))
                  end do 

               end do

            end if

         end do

      end do

      end

      subroutine makapc (id)
c----------------------------------------------------------------------
c subroutine to make the ap'*c matrix for ap'*c = c, where p' is the 
c the independent subset of p.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, m, n, id

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k10)

      double precision apc
      common/ cstp2c /apc(h9,k5,m19)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer jend
      common/ cxt23 /jend(h9,m4)
c----------------------------------------------------------------------
      apc(id,:,:) = 0d0
c                                 independent endmembers:
      do j = 1, lstot(id)
c                                 the disordered endmembers
         do i = 1, icomp
            apc(id,i,j) = cp(i,jend(id,2+j))
         end do

      end do

      if (lstot(id).lt.nstot(id)) then 
c                                  if o/d add the ordered endmembers
         do m = 1, norder
c                                  the ordered endmember compositions are 
c                                  a stoichiometric combination of the 
c                                  disordered endmembers
            j = lstot(id) + m

            do n = 1, nr(m)
c                                  the endmember decomposes to iddeps(n,m) in the
c                                  y array, or ideps(n,m,id) in the p array.
               do i = 1, icomp

                  apc(id,i,j) = apc(id,i,j) + depnu(n,m)
     *                        * cp(i,jend(id,2+ideps(n,m,id)))
               end do 

            end do

         end do

      end if
c                                   eliminate p(nstot) as 1 - sum(p,nstot-1)
c                                   this method costs an extra column dimension in
c                                   p2c, but what the heck... it can be used 
c                                   to form the constraint b vector
      do j = 1, nstot(id) - 1
         do i = 1, icomp
            apc(id,i,j) = apc(id,i,j) - apc(id,i,nstot(id))
         end do
      end do

      end

      subroutine setord (im)
c---------------------------------------------------------------------
c set global order/disorder models parameters, call by gmodel for 
c solution model im. 
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical zbad

      integer im, i, j, ind, id, k, l,itic, ii, imatch, 
     *        il, ik, kk, jp1

      double precision dzt, dx, delta, c0(0:20), c1(0:20), zij

      external zbad

      character tname*10
      logical refine, resub
      common/ cxt26 /refine,resub,tname

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision p,t,xco2,mmu,tr,pr,r,ps
      common/ cst5 /p,t,xco2,mmu(2),tr,pr,r,ps

      integer nsub,nterm
      double precision acoef
      common/ cst107 /acoef(m10,m11,0:m0),
     *                nterm(m10,m11),nsub(m10,m11,m0)

      integer iend,isub,insp,iterm,iord,istot,jstot,kstot,rkord
      double precision wg,wk
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),iend(m4),
     *      isub(m1,m2),insp(m4),
     *      rkord(m18),iterm,iord,istot,jstot,kstot

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision deph,dydy,dnu
      common/ cxt3r /deph(3,j3,h9),dydy(m4,j3,h9),dnu(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m4)

      double precision cp
      common/ cst12 /cp(k5,k10)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)

      character specie*4
      integer jsp, ins
      common/ cxt33 /jsp,ins(nsp),specie(nsp)

      character mname*8
      common/ cst18a /mname(m4)

      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
c                                 models with speciation:
      do j = 1, norder

         do i = 1, 3
            deph(i,j,im) = denth(j,i)
         end do

         nrct(j,im) = nr(j)

         do i = 1, nr(j)
            ideps(i,j,im) = iy2p(iddeps(i,j))
         end do

      end do
c                                 classify multiple species models according
c                                 to whether the disordered reactants are
c                                 partially or completely correlated, assume
c                                 anti-correlation is not a possible case.
      icase(im) = 0

      if (norder.gt.1) then

         imatch = 0

         do j = 1, nr(1)

            id = ideps(j,1,im)

            do i = 1, nr(2)
               if (id.eq.ideps(i,2,im)) then
                  imatch = imatch + 1
                  exit
               end if
            end do

         end do

         if (imatch.eq.1) then
c                                 if match = 1 one species didn't match
c                                 assume partial correlation
            icase(im) = 2
         else if (imatch.ge.2) then
            icase(im) = 1
         end if

      end if
c                                first create derivatives of endmember
c                                fractions with respect to the ordered
c                                species:
      do j = 1, norder

         do i = 1, nstot(im)
            dydy(i,j,im) = 0d0
         end do
c                                derivative of the ordered species with
c                                respect to itself:
         dydy(kstot+j,j,im) = 1d0
c                                each ordered species decomposes to
c                                two disordered species iddeps(1-2,j)
c                                depnu is the stoichiometric coefficient
c                                of the disordered species in the ordered
c                                species.

c                                derivatives of the consituent species
c                                with respect to the ordered species
         dnu(im) = 1d0

         do i = 1, nr(j)
            dydy(ideps(i,j,im),j,im) = dydy(ideps(i,j,im),j,im)
     *                                  - depnu(i,j)
            dnu(im) = dnu(im) + dydy(ideps(i,j,im),j,im)
         end do
c                                dnu ~0 => speciation reaction is not equimolar
         if (dabs(dnu(im)).gt.1d4*zero) then
            if (norder.gt.1) call error (72,r,i,
     *              'ordering schemes with > 1 non-equi'//
     *              'molar reaction have not been anticipated: '//tname)
         else
            dnu(im) = 0d0
         end if

      end do
c                                evaluate the second derivative of each
c                                pi*pj term in the excess function with
c                                respect to kth species
      do i = 1, iterm
         do j = 1, norder
            do k = 1, norder

                  dppp(k,j,i,im) =  dydy(jsub(1,i,im),k,im)
     *                             *dydy(jsub(2,i,im),j,im)
     *                           +
     *                              dydy(jsub(2,i,im),k,im)
     *                             *dydy(jsub(1,i,im),j,im)
            end do
         end do
      end do
c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      do i = 1, msite(im)
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         if (zsp(im,i)+1.gt.m11) call error (1,dx,zsp(im,i)+1,'m11')

         do k = 1, norder
            sdzdp(k,zsp(im,i)+1,i,im) = 0d0
         end do

         do j = 1, zsp(im,i)
c                                 # of terms in the
c                                 site fraction function and a0.
            do l = 1, norder
               sdzdp(l,j,i,im) = 0d0
            end do
c                                 for each term:
            do k = 1, lterm(j,i,im)
c                                 endmember indexes
               ind = ksub(k,j,i,im)
c                                 get derivatives, of species fractions
c                                 with respect to ordered species
               do l = 1, norder
                  itic = 0
                  do ii = 1, nr(l)
                     if (ind.eq.ideps(ii,l,im)) then
                        sdzdp(l,j,i,im) = sdzdp(l,j,i,im)
     *                  + dydy(ideps(ii,l,im),l,im)*dcoef(k,j,i,im)
                        itic = itic + 1
c                                 high order terms not allowed
                        if (itic.gt.1) call error (999,r,801,'GMODEL')
                     end if
                  end do
c                                 the derivative of a term with the
c                                 ordered species.
                  if (ind.eq.kstot+l)
     *               sdzdp(l,j,i,im) = sdzdp(l,j,i,im)
     *                               + dcoef(k,j,i,im)

               end do
            end do
         end do
      end do
c                                 multiply each dzdp by qmult (R*site
c                                 multiplicity) to reduce operation
c                                 count in evaluating derivatives.
      do k = 1, norder
         do i = 1, msite(h0)

            dzt = 0d0

            do j = 1, zsp(im,i)
               if (dabs(sdzdp(k,j,i,im)).lt.zero) sdzdp(k,j,i,im) = 0d0
               dzt = dzt + sdzdp(k,j,i,im)
            end do

            if (dabs(dzt).lt.zero) dzt = 0d0
            sdzdp(k,j,i,im) = -dzt

         end do
      end do
c                                 ----------------------------------------------
c                                 derive z2p limit expressions for O/D models
      do k = 1, nord(im)
c                                 number of limits for ordered species k
         ln(k,im) = 0
      end do

      if (ksmod(im).ne.688) then 

         if (nord(im).gt.1) call error (72,c0(0),i,
     *            'solution '//tname//': multiple order parameters '//
     *            'only allowed in 688 format models')

c                                 all z expressions may be necessary to
c                                 formulate limits, make the ksp'th + 1
c                                 species expression by differnce
         do i = 1, msite(im)
c                                 qmult = 0, temkin, all expressions are
c                                 available
            if (zmult(im,i).eq.0) cycle

            jp1 = zsp(im,i) + 1
c                                 initialize the term counter
            lterm(jp1,i,im) = 0
c                                 cycle through the endmembers to work out
c                                 if it has a non zero fraction
            do l = 1, nstot(im)

               call zmake (zij,i,l,im)

               if (dabs(zij).lt.zero) cycle

               lterm(jp1,i,im) = lterm(jp1,i,im) + 1
               dcoef(lterm(jp1,i,im),jp1,i,im) = zij
               ksub(lterm(jp1,i,im),jp1,i,im) = l

            end do
         end do
      end if 

      do i = 1, msite(im)

         if (zmult(im,i).eq.0) then
            jp1 = 0
         else
            jp1 = 1
         end if

         do j = 1, zsp(im,i) + jp1

            do k = 1, nord(im)

               c0 = 0d0
               c1 = 0d0

               if (dcoef(0,j,i,im).ne.0d0) call error (72,c0(0),i,
     *            'solution '//tname//': constants not allowed in '//
     *            'O/D model configurational entropy site fraction '//
     *            'expressions')

               do l = 1, lterm(j,i,im)

                  il = ksub(l,j,i,im)

                  if (il.le.lstot(im)) then

                     c0(il) = c0(il) + dcoef(l,j,i,im)

                     do ik = 1, nord(im)

                        kk = lstot(im) + ik
c                                  coefficient on p0
                        c0(kk) = c0(kk) - dydy(il,ik,im)
     *                                  * dcoef(l,j,i,im)
c                                  coefficient on p
                        c1(kk) = c1(kk) + dydy(il,ik,im)
     *                                  * dcoef(l,j,i,im)

                     end do

                  else

                     c1(il) = c1(il) + dcoef(l,j,i,im)

                  end if
c                                 at this point the c's are the coefficients for a
c                                 z(p,p0), below they are rearranged to get p(kk) = f(p0,p[~kk],z[0,1])
c                                 in other words the loop on mord(im) is superfluous, but what the heck...
                  kk = lstot(im)+k

                  if (l.eq.lterm(j,i,im).and.c1(kk).ne.0d0) then

                     do ik = 0, nstot(im)
                        c0(ik) = -c0(ik)/c1(kk)
                     end do
c                                the constant for the p(k) limit when z(j) = 1
                     c1(0) = c0(0) + 1d0/c1(kk)

                     do ik = lstot(im) + 1, nstot(im)
                        if (ik.eq.kk) cycle
                        c1(ik) = -c1(ik)/c1(kk)
                     end do

                     c1(kk) = 0d0

                     if  (c1(0).gt.c0(0)) then
c                                 z = 1 is the upper limit, the constant is c0(0) and
                        delta = c1(0)-c0(0)

                     else
c                                 z = 0 is the upper limit, the constant is c1(0) and
                        delta = c0(0)-c1(0)
                        c0(0) = c1(0)

                     end if
c                                 -----------------------------------------------------
c                                 found a limit, increment limit counter
                     ln(k,im) = ln(k,im) + 1
c                                 initialize p0 term counter for limit
                     lt(ln(k,im),k,im) = 0
c                                 load the p0 coefficients, if simplicial p0 > lstot(im) = 0
                     do ik = 1, nstot(im)

                        if (jsmod.eq.6.and.ik.gt.lstot(im)) exit
                        if (c0(ik).eq.0d0) cycle
c                                 increment term counter:
                        lt(ln(k,im),k,im) = lt(ln(k,im),k,im) + 1
c                                 save the coefficient and index:
                        lc(lt(ln(k,im),k,im),ln(k,im),k,im) = c0(ik)
                        lid(lt(ln(k,im),k,im),ln(k,im),k,im) = ik

c                          write (*,*) i,j,ln(k,im),lt(ln(k,im),k,im)
c                          write (*,*) c0(ik),ik

                     end do
c                                 initialize p term counter for limit
                     jt(ln(k,im),k,im) = 0
c                                 load the p coefficients
                     do ik = lstot(im) + 1, nstot(im)

                        if (c0(ik).eq.0d0) cycle
c                                 increment term counter:
                        jt(ln(k,im),k,im) = jt(ln(k,im),k,im) + 1
c                                 save the coefficient and index:
                        jc(jt(ln(k,im),k,im),ln(k,im),k,im) = c1(ik)
                        jid(jt(ln(k,im),k,im),ln(k,im),k,im) = ik

c                          write (*,*) i,j,jt(ln(k,im),k,im)
c                          write (*,*) c0(ik),ik


                     end do
c                                load the constant and delta:
                     l0c(1,ln(k,im),k,im) = c0(0)
                     l0c(2,ln(k,im),k,im) = delta

c                       write (*,*) 'cst delta ', c0(0),delta
c                       write (*,*) ' '

                  end if

               end do

            end do

         end do

      end do

      end

      subroutine minfxc (ids,maxs)
c-----------------------------------------------------------------------
c optimize solution gibbs energy or configurational entropy at constant 
c composition subject to site fraction constraints.

c     number of compositional constraints -> icomp (<k5)
c     number of independent endmember fractions -> nstot-1 (<m19)
c     number of site fractions -> nz (<m20)
c     closure is forced in the objective function (gsol2) and
c        by using the complete set of site fraction.
c     requires that pp has been loaded in cxt7
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical maxs

      integer ids, i, j, k, nvar, iter, iwork(m22), iprint,
     *        istuff(10),istate(m21), istart, nclin, ntot

      double precision sum, ggrd(m19), b(k5),
     *                 bl(m21), bu(m21), gfinal, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),stuff(1),
     *                 lapz(m20,m19)

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision apc
      common/ cstp2c /apc(h9,k5,m19)

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

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      external gsol3
c-----------------------------------------------------------------------
c                                 create constraint matrix z
c                                 --------------------------
      ntot = nstot(ids)
      nvar = ntot - 1
c                                 initialize bounds
      bu(1:nvar) = 1d20
      bl(1:nvar) = -1d20
c                                 load the local constraints 
c                                 from the global arrays
      nclin = nz(ids) + icomp
c                                 first the site fraction constraints
      lapz(1:nclin,1:nvar) = apz(ids,1:nclin,1:nvar)

      do i = 1, nclin
         bl(nvar+i) = zl(ids,i)
         bu(nvar+i) = zu(ids,i)
      end do
c                                 get the bulk composition from pp
      call getscp (b,sum,ids,ids,.false.)
c                                 the bulk composition constraints
      do k = 1, icomp

         i = nz(ids) + k

         do j = 1, nvar
            lapz(i,j) = apc(ids,k,j)
         end do

         bl(nvar+i) = b(k) - apc(ids,k,ntot)
         if (dabs(bl(nvar+i)).lt.zero) bl(nvar+i) = 0d0
         bu(nvar+i) = bl(nvar+i)

      end do

      istart = -1
c                                 solution model index
      istuff(1) = ids
c                                 istuff(2) is set by NLP and is
c                                 irrelevant.
c                                 istuff(3) = 0 min g, 1 max entropy
      if (maxs) then 
         istuff(3) = 0
      else 
         istuff(3) = 1
      end if

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

c                                 istuff = 0 min g, 1 max entropy
      if (maxs) then 
         istuff(1) = 0
      else 
         istuff(1) = 1
      end if 

      call e04ucf (nvar,nclin,m20,m19,lapz,bl,bu,gsol3,
     *             iter,istate,clamda,gfinal,ggrd,r,ppp,iwork,
     *             m22,work,m23,istuff,stuff,istart,iprint)

      end
