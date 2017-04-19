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
     *        /,'Ionic strength = ',g12.6,' gamma/q^2 = ',g12.6,
     *        'Permativity =',g12.6,//,10x,'  molality ')
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