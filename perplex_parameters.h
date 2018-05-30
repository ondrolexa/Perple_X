

      integer h5,h6,h8,h9
      integer i6,i7,i8,i9,i10,i11
      integer j3,j4,j5,j6,j9
      integer k0,k1,k2,k3,k4,k5,k7,k8,k9,k10,k13,k14,k15
      integer k16,k17,k18,k19,k20,k21,k22,k23,k24,kd2
      integer l2,l3,l5,l6,l7,l8,l9,l10,lchar
      integer m0,m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15
      integer m16,m17,m18
      integer msp,mst,mdim,ms1
      integer n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,nsp,nx,ny
!                                 n0  - starting LUN-1 for fractionation files, these files may 
!                                       have LUNs from n0+1 up to n0+k23
!                                 n1  - problem definition file.
!                                 n2  - thermodynamic data file.
!                                 n3  - print output file
!                                 n4  - plot output file 1
!                                 n5  - plot output file 2 (bplot)
!                                 n6  - reaction list file
!                                 n7  - unused
!                                 n8  - locally opened and closed files
!                                 n9  - solution model file 
!                                 n10 - autorefine file 1
!                                 n11 - autorefine file 2
!                                 n12 - failed result file (fld)
!                                 n13 - aq error file (pts)
      parameter (n10=7,n11=8,n1=11,n2=12,n3=13,n4=14,n5=15,n6=16,n7=17)
      parameter (n8=18,n9=19,n12=20,n13=21,n0=30)
!                                 msp - max number of species on a solution identisite
!                                 mst - max number of distinct identisites per solution
!                                 mdim - hard constraint on max number of dimensions
!                                        for a solution model composition space.
      parameter (mst=3,mdim=8,msp=mdim+2,ms1=msp-1)
!                                 h5 - max number of saturated components
!                                 h6  - max number of saturated composants in any subcomposition
!                                 h8  - max number of excluded phases
!                                 h9  - max number of solutions
      parameter (h5=5,h6=500,h8=200,h9=30)
!                                 i6  - maximum number of independent chemical potentials (or 
!                                       fugacity/activities).
!                                 i7  - number of system props used in werami
!                                 i8  - number of properties saved for each phase and
!                                       for the bulk composite, less than i8 properties
!                                       may be read from the bulk property file. 
!                                       if i8 is changed then adjust:
!                                          pt2prp - frendly.f
!                                          l2p    - werami.f
!                                          prname - werami.f
!                                 i9  - max number of solutions in solution model file
!                                 i10 - max number of option values in perplex_option.dat
!                                 i11 - max number of dependent properties in a tab format file
      parameter (i6=2,i7=20,i8=28,i9=200,i10=50,i11=100)
!                                 j3  - max number of ordered species
!                                 j4  - max number of species in the definition of a dependent species
!                                 j5  - max number of stoichiometric limits on an ordered species
!                                 j6  - max number of terms in a stoichiometric limit on an ordered species
!                                 j9  - max number of divariant assemblages
      parameter (j3=4,j4=6,j5=8,j6=8,j9=160000)
!                                 k0  - max number of database components
!                                 k1  - max number of compounds
!                                 k2  - max number of invariant and univariant compound 
!                                       equilibria and max number of divariant compound
!                                       assemblages for constrained bulk composition
!                                       calculations
!                                 k3  - max number of distinct phase (as opposed to 
!                                       pseudocompound assemblages)
!                                 k4  - number of parameters in a Perple_X EoS
!                                 k5  - max total number of active components
!                                 k9  - maximum number of true compounds with lambda transitions
!                                 k10 - max number of true compounds
!                                 k13 - pseudocompound parameter array dimension
!                                 k14 - number of parameters in the data base EoS
!                                 k15 - number of elastic moduli parameters in the emod array
!                                 k16 - max number of make definitions
!                                 k17 - max number of entities in a make definition
!                                 k18 - max number of "x" coordinates saved for pseudocompunds,
!                                       this cannot exceed (k1-k10)*mst*(msp-1), but in general
!                                       will be much smaller (this parameter is similar to k13).
!                                 k19 - max number of refinement points for adaptive optimization.
!                                 k20 - max number of "z" coordinates saved for pseudocompounds
!                                       generated by adaptive refinement, max value = k21*m4, 
!                                       but usually much smaller (see k18).
!                                 k21 - max number of pseudocompounds for adaptive refinement.
!                                 k23 - max number of phases to be fractionated.
!                                 k24 - max number of coordinates for the pseudocompounds of a 
!                                       prismatic solution.
      parameter (k0=25,k1=1500000,k2=100000,k3=2000,k4=32,k5=12)
      parameter (k7=k5+1,k8=k5+2) 
      parameter (k9=30,k10=400,k14=18,k15=6,k16=100)
      parameter (k17=7,k18=29*k1)
      parameter (k19=3*k5,k21=2000000,k20=(mdim+3)*k21)
      parameter (k22=mdim*k19,k23=25)
      parameter (k13=mdim*k21,k24=k13*(mst-1))
!                                 l2 - max number of independent potential variables
!                                 l3 - max number of variables for gridded min and graphics (l2+2)
!                                 l5 - max number of coordinates along a univariant curve                
!                                 l6 - max number of iterations in lp optimization, in 
!                                      theory this may be up to ca 5*(k1+k5), generally
!                                      convergence occurs in less than 100 iterations.
!                                 l7 - max number of grid points along an axis for
!                                      constrained minimization on a 2-d grid. this
!                                      array is not essential.  
!                                 l8 - max number of levels for multilevel grids    
!                                 l9 - max number of aqueous solute species in minimization programs.         
!                                l10 - max number of parameters stored in caq for each phase.
!                                nsp - max number of species in fluid speciation routines 
      parameter (l2=5,l3=l2+2,l5=1000,l6=1000,l7=2048,l8=10,l9=100,
     *           nsp=17,l10=nsp+l9+4)
!                                 m0 - max number of terms for a species site fraction?
!                                 m1 - max number of terms in excess function
!                                 m2 - max order of term in excess function
!                                 m3 - number of parameters for excess function coefficients
!                                 m4 - max number of endmembers in a solution model,
!                                      should include dependent endmembers? this cannot
!                                      exceed msp**mst+1, but will generally be MUCH
!                                      smaller.
!                                 m6 - max number of "lambda" transitions per true compound
!                                 m7 - number of "lambda" parameters per transition
!                                 m8 - maximum number of parameters to describe T dependent ordering.
!                                 m9 - maximum number of true compounds with T-dependent ordering.
!                                m10 - maximum number of sites in Sconf model
!                                m11 - maximum number of species - 1 in Sconf model
!                                m12 - maximum order of terms in Sconf model
!                                m13 - maximum number of terms in endmember fraction
!                                      function expressed as a function of dependent
!                                      endmember
!                                m14 - max order of term in endmember fraction
!                                      function expressed as a function of dependent
!                                      endmember
!                                m15 - maximum number of dependent endmembers for a 
!                                      reciprocal solution
!                                m16 - max number of parameters in a redkich-kistler L
!                                m17 - max order of redlich-kistler expansion
!                                m18 - max number of pairwise terms in a redlich-kistler expansion
      parameter (m0=8,m1=60,m2=8,m3=3,m4=20,m6=6,m7=15,m8=9,m9=10,m10=6,
     *           m11=6,m12=4,m13=8,m14=2,m15=12,m16=6,m17=5,m18=6)
!                                 nx - number of x-grid nodes in a contour data grid
!                                 ny - number of y-grid modes in a contour data grid
      parameter (nx=500,ny=500)
c                                 lchar - maximum length of character strings
      parameter (lchar=400)

! NOTE: increasing parameter K5 requires changes to the following
! format statements:
!                     1000 in routine bulkck
!                     1000 in routine bulktst
!                     1150 in routine input1
!                     4000 in main in program build


!                                  the following are dependent parameters (also k7,k8)
      parameter (kd2=k8*35)

