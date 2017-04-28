 
      subroutine solvent (ins,isp)
c-----------------------------------------------------------------------
c compute the gibbs energies, densities and dielectric constants for the
c pure molecular species of a solvent
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j

      double precision po(nsp,8), eps(k)
c                                Harvey & Lemmon provide additional data for 
c                                ethylene and long-chain hydrocarbons. A_mu
c                                is zero for all species listed here, therefore
c                                Eq 5 of H&L:
c                                P/rho (cm3/mol) = A + A_mu/T + B*rho + C*rho^D
c                                is simplified by dropping the second term the 
c                                coefficients are a f(T) viz A = a0 + a1*(T/Tr - 1)....
c                             
c                                por(i,1:8) - a0, a1, A_mu, b0, b1, c0, c1, D
      data por(i,j),j=1,8),i=1,nsp)/ 
c                                1 - H2O
     *     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 
c                                2 - CO2
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0, 
     *     1.55d0,
c                                3 - CO approximated by O2
     *     3.9578d0, 6.5d-3, 0d0, 0.575d0, 1.028d0, -8.96d0, -5.15d0, 
     *     1.5d0,
c                                4 - CH4
     *     6.5443d0, 1.33d-2, 0d0,8.4578d0, 3.7196d0, -352.97d0, 
     *     -100.65d0, 2d0,
c                                5 - H2
     *     2.0306d0, 5.6d-3, 0d0, 0.181d0, 0.021d0, -7.4d0, 0d0, 2d0,
c                                6 - H2S approximate by ?
     *     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 
c                                7 - O2
     *     3.9578d0, 6.5d-3, 0d0, 0.575d0, 1.028d0, -8.96d0, -5.15d0, 
     *     1.5d0,
c                                8 - SO2 approximated by CO2
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0, 
     *     1.55d0,
c                                9 - COS approximated by CO2
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0, 
     *     1.55d0,
c                                10 - N2
     *     4.3872d0, 2.26d-3, 0d0, 2.206d0, 1.135d0, -169d0, -35.83d0, 
     *     2.1d0,
c                                11 - NH3 approximated by ?
     *     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
c                                12-15 Si-O high T species
     *     32*0d0
c                                16 - Ethane
     *     11.1552d0, 0.0112d0, 0d0, 36.759d0, 23.639d0, -808.03d0, 
     *     -378.84d0, 1.75/

      save por 
c----------------------------------------------------------------------
      rhoi = 1d0
      trt = t/tr - 1d0

      do i = 1, ns - 1

         j = ins(i)
c                                 Eq 5 of H&L 2005 for pj = polarization/rho
c                                 for polar species need to add
         pj = po(j,1) + po(j,2)*trt + (po(j,4) + po(j,5)*trt)*rhoi
        *             + (po(j,6) + po(j,7)*trt)*rhoi**po(j,8)
c                                 invert clausius-mosotti relation for dielectric
c                                 constant (Eq 1) non-polar molecules
         eps(i) = (rhoi - 2d0*pj) / (pj - rhoi)
c                                 invert kirkwood relation for dielectric
c                                 constant (Eq 2), this would be necessary
c                                 if H2S were present.            
c        eps(i) = (9d0*pj + rhoi 
c     *          + 3d0*dsqrt(9d0*pj*pj+2d0*pj*rhoi+rhoi*rhoi)) / (9d0*pj)

      end do 

      end 


      subroutine aqrxdo
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, l, ichg, ihy, jchg(l9), it

      logical bad, output

      double precision c(l9), q(l9), mo(l9), dg(l9), gh2o, ahy, xis,
     *                 d(l9), q2(l9), lng0, is, gamm0

      double precision gcpd, solve

      external gcpd, solve

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision vh2o, epsilo, adh
      common/ cxt37 /vh2o, epsilo, adh

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      output = .true. 
c                                 free energies of aqueous species
      ichg = 0 
      is = 0d0
c DEBUG DEBUG
      gh2o = mu(1)

      do i = 1, aqct 

         k = aqst + i 

         if (thermo(1,k).eq.0d0) then 
c                                  find hyrdonium
            ihy = i 
            cycle 

         end if 

         q(i) = thermo(6,k)
         q2(i) = q(i)*q(i)
c                                 dg is the solvent oxide potentials - g
         dg(i) = -gcpd(k,.false.)

         do j = 1, jbulk 
c                                 if oxide components, but no excess oxygen
c                                 mu(O2) is a nan. 
            if (isnan(mu(j))) cycle

            dg(i) = dg(i) + aqcp(j,i) * mu(j)

         end do 
c                                 normalize by RT
         dg(i) = (dg(i) - q(i)*gh2o/2d0)/r/t

         if (q(i).ne.0d0) then 

            ichg = ichg + 1
            jchg(ichg) = i
c                                  gh2o is the partial molar gibbs energy of
c                                  water, don't use mu because we don't know 
c                                  if h2o is a component. 
c                                  this is now c(i)*a(H+)^(q(i)) = mol(i)*gamma(i)*q(i)
            d(i) = q(i)*dexp(dg(i))
            c(i) = d(i)

         else 
c                                  neutral species assumed to be ideal, molality is
            mo(i) = dexp(dg(i))
      
         end if 

      end do 

      gamm0 = 1d0 
      it = 0 
      mo(ihy) = 1e-6
c                                  iterative loop for ionic strength
      do 
c                                  solve charge balance for H+
         mo(ihy) = solve(c,q,mo(ihy),jchg,ichg,bad)

         if (bad) then 
            write (*,*) 'bombed'
            pause
         end if 

         ahy = mo(ihy)*gamm0
c                                  back calculate charged species molalities
c                                  and ionic strength
         xis = is

         is = 0d0

         do i = 1, ichg 

            j = jchg(i)
            mo(j) = ahy**(q(j))*c(j)/q(j)
            is = is + q2(j) * mo(j)

         end do

         is = is / 2d0 
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
         lng0 = adh*dsqrt(is)/(1d0 + dsqrt(is)) + 0.2d0*is
         gamm0 = dexp(lng0)
c                                 check for convergence
         if (dabs(xis-is)/is.lt.nopt(5)) then 
            exit
         else if (it.gt.iopt(21)) then 
            bad = .true.
         end if 

         it = it + 1
c                                 update coefficients
         do i = 1, ichg 

            j = jchg(i)
            c(j) = d(j)*dexp(lng0*(1d0-q2(j)))

         end do

      end do

      if (output) then
         write (*,1000) is,gamm0,epsilo
         do i = 1, aqct
            write (*,1010) aqnam(i),mo(i)
         end do 
      end if

1000  format (/,'Rock-dominated solvent solute speciation:',/,
     *        /,'Ionic strength = ',g12.6,'; gamma/q^2 = ',g12.6,
     *        '; Permativity =',g12.6,//,10x,'  molality ')
1010  format (a8,2x,g12.6)

      end 

      double precision function solve (c,q,x,jchg,ichg,bad)
c-----------------------------------------------------------------------
c function solve for hydronium molality (x) from charge balance for aqrxdo.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer jchg(*), ichg, it, i, j

      logical bad

      double precision c(*), q(*), x, y, z, f, df, dx

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------

      it = 0 

      do

         f = x
         df = 1d0
         it = it + 1

         do i = 1, ichg 

            j = jchg(i)

            y = x**(q(j)) * c(j)
            z = y*q(j)/x 

            f = f + y
            df = df + z

         end do  

         dx = -f/df

         x = x + dx

         if (dx/x.lt.nopt(5)) then 
            bad = .false.
            exit
         else if (x.lt.0d0.or.it.gt.iopt(21)) then 
            bad = .true.
            exit
         end if 

      end do

      solve = x

      end 

      subroutine setp0a (ids,id)
c-----------------------------------------------------------------------
c for static pseudocompounds load the compositional coordinates from xco
c into simple compositional arrays for compound id of solution ids. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ids, id, i

      double precision xt
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /x(m4),y(m4),pa(m4),p0a(m4),z(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)
c-----------------------------------------------------------------------

      xt = 0d0 

      if (lorder(ids).and.lrecip(ids)) then 

         do i = 1, nstot(ids) - 1
            p0a(i) = xco(jco(id)+i) 
            pa(i) = p0a(i)
            xt = xt + pa(i)
         end do 

         p0a(i) = 1d0 - xt  
         pa(i) = p0a(i)
  
      else if (lorder(ids)) then

         do i = 1, lstot(ids) - 1
            p0a(i) = xco(jco(id)+i) 
            pa(i) = p0a(i)
            xt = xt + pa(i)
         end do

         p0a(i) = 1d0 - xt
         pa(i) = p0a(i)

      else 

         do i = 1, lstot(ids)- 1
            p0a(i) = xco(jco(id)+i) 
            xt = xt + p0a(i)
         end do

         p0a(i) = 1d0 - xt

       end if

       end