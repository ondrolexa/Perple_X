 
      subroutine slvnt1
c-----------------------------------------------------------------------
c computes dielectric constant for the molecular species of a solvent
c assumes species volumes or partial molar volumes in cohhyb have been 
c initialized.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j

      double precision po(nsp,8), pj, trt, rho

      double precision epsh2o

      external epsh2o

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision y,g,v,eps,v0,eps0
      common/ cstcoh /y(nsp),g(nsp),v(nsp),eps(nsp),v0(nsp),eps0(nsp)
c                                Harvey & Lemmon provide additional data for 
c                                ethylene and long-chain hydrocarbons. A_mu
c                                is zero for all species listed here, therefore
c                                Eq 5 of H&L:

c                                P/rho (cm3/mol) = A + A_mu/T + B*rho + C*rho^D

c                                is simplified by dropping the second term the 
c                                coefficients are a f(T) viz A = a0 + a1*(T/Tr - 1)...

c                                po(i,1:8) - a0, a1, A_mu, b0, b1, c0, c1, D
      data ((po(i,j),j=1,8),i=1,nsp)/ 
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
     *     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,
c                                12-15 Si-O high T species
     *     32*0d0,
c                                16 - Ethane
     *     11.1552d0, 0.0112d0, 0d0, 36.759d0, 23.639d0, -808.03d0, 
     *     -378.84d0, 1.75d0,
c                                17 - dilutant
     *      8*0d0/

      save po 
c----------------------------------------------------------------------

      trt = t/273.16d0 - 1d0

      do i = 1, ns - 1

         j = ins(i)
c                                 rho = 1/vcm3, v(j) is initialized by lnfpur
c                                 in gcpd and is in j/bar/mol
         rho = 0.1d0/v(j)
c                                 Eq 5 of H&L 2005 for pj = polarization/rho
c                                 for polar species need to add
         pj = po(j,1) + po(j,2)*trt + (po(j,4) + po(j,5)*trt)*rho
     *                + (po(j,6) + po(j,7)*trt)*rho**po(j,8)

         if (po(j,3).eq.0d0) then 
c                                 invert clausius-mosotti relation for dielectric
c                                 constant (Eq 1) non-polar molecules
            eps(i) = (2d0*pj*rho + 1d0) / (1d0 - rho*pj)
         else
c                                 polarized species, currently none, but H2S...
            pj = pj + po(j,3)

            eps(i) = (9d0*pj*rho + 1d0 
     *             + 3d0*dsqrt((3d0*pj*rho)**2+2*p*rho+1))/4d0

         end if

      end do 
c                                 and water:
      eps(i) = epsh2o (v(ins(i))) 

      end 

gmodel:

         if (ksmod(im).eq.20) then 
c                                 compute endmember reference permittivities 
c                                 for HKF electrolyte model
           p = pr
           t = tr
c                                 this is just to set the volumes
           do i = 1, ns - 1
c              dzt = gcpd(jnd(i),.true.)
              v(ins(i)) = r*t/p
           end do 

           v(ins(ns)) = 1.80685d0

           call slvnt1

           do i = 1, ns
              eps0(i) = eps(i)
              v0(ins(i)) = v(ins(i))
           end do 

         end if 