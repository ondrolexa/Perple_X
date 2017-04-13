      subroutine aqrxdo
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, l, ichg, ihy, jchg(l9)

      double precision c(l9), q(l9), mo(l9), dg(l9), gh2o, f, df, x, 
     *                 dx, y, z

      double precision gcpd

      external gcpd 

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer nq,nn,ns,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,nqs,nqs1,sn,qn,nq1

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
c-----------------------------------------------------------------------
c                                 free energies of aqueous species
      ichg = 0 
c DEBUG DEBUG
      gh2o = mu(1)/r/t

      do i = 1, aqct 

         k = aqst + i 

         if (thermo(1,k).eq.0d0) then 
c                                  find hyrdonium
            ihy = i 
            cycle 

         end if 

         q(i) = thermo(6,k)
c                                 dg is the solvent oxide potentials - g
         dg(i) = -gcpd(k,.false.)/r/t

         do j = 1, jbulk 
c                                 if oxide components, but no excess oxygen
c                                 mu(O2) is a nan. 
            if (isnan(mu(j))) cycle

            dg(i) = dg(i) + aqcp(j,i) * mu(j)

         end do 
c                                 normalize by RT
         dg(i) = dg(i)/r/t

         if (q(i).ne.0d0) then 

            ichg = ichg + 1
            jchg(ichg) = i
c                                  gh2o is the partial molar gibbs energy of
c                                  water, don't use mu because we don't know 
c                                  if h2o is a component. 
c                                  this is now c(i)*a(H+)^(q(i)) = mol(i)*gamma(i)*q(i)
            c(i) = q(i)*dexp(dg(i)- q(i)*gh2o/2d0)

         else 
c                                  neutral species assumed to be ideal, molality is
            mo(i) = dexp(dg(i))
      
         end if 

      end do 
c                                  solve charge balance for H+
      x = 1e-4

      do

         f = x
         df = 1d0 

         do i = 1, ichg 

            j = jchg(i)

            y = x**(q(j)) * c(j)
            z = y*q(j)/x 

            f = f + y
            df = df + z

         end do  

         dx = -f/df

         x = x + dx
         if (x.lt.0d0) then 
            write (*,*) 'oink'
         end if 

      end do

      end 
