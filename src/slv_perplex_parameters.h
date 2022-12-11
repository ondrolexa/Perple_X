! special dimensions for SLV calcs
      integer h5,h6,h8,h9
      integer i6,i7,i8,i9,i10,i11
      integer j3,j4,j5,j6,j9
      integer k0,k1,k2,k3,k4,k5,k7,k8,k9,k10,k12,k13,k14,k15
      integer k16,k17,k18,k19,k20,k21,k22,k23,k24,kd2
      integer l2,l3,l5,l6,l7,l8
      integer m0,m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15
      integer msp,mst,nsp,mdim,ms1
      integer n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,nx,ny
!                                 msp - max number of species on a solution identisite
!                                 mst - max number of distinct identisites per solution
!                                 mdim - hard constraint on max number of dimensions
!                                        for a solution model composition space.
      parameter (msp=9,mst=2,mdim=8,ms1=msp-1)
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
      parameter (n10=7,n11=8,n1=11,n2=12,n3=13,n4=14,n5=15,n6=16,n7=17)
      parameter (n8=18,n9=19,n0=30)
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
!                                 i9  - max number of solution in solution model file
!                                 i10 - max number of option values in perplex_option.dat
!                                 i11 - max number of dependent properties in a tab format file
      parameter (i6=2,i7=20,i8=27,i9=200,i10=25,i11=10)
!                                 j3  - max number of ordered species
!                                 j4  - max number of species in the definition of a dependent species
!                                 j5  - max number of stoichiometric limits on an ordered species
!                                 j6  - max number of terms in a stoichiometric limit on an ordered species
!                                 j9  - max number of divariant assemblages
      parameter (j3=4,j4=5,j5=6,j6=8,j9=1)
!                                 k0  - max number of database components
!                                 k1  - max number of compounds
!                                 k2  - max number of invariant and univariant compound 
!                                       equilibria and max number of divariant compound
!                                       assemblages for constrained bulk composition
!                                       calculations
!                                 k3  - max number of distinct phase (as opposed to 
!                                       pseudocompound assemblages)
!                                 k4  - number of parameters in PerpleX EoS
!                                 k5  - max total number of active components
!                                 k9  - maximum number of true compounds with lambda transitions
!                                 k10 - max number of true compounds
!                                 k12 - max number of end-members in a solution model
!                                 k13 - pseudocompound parameter array dimension
!                                 k14 - number of parameters in the data base EoS
!                                 k15 - number of elastic moduli parameters in the emod array
!                                 k16 - max number of make definitions
!                                 k17 - max number of entities in a make definition
!                                 k18 - max number of "x" coordinates saved for pseudocompunds,
!                                       this cannot exceed (k1-k10)*mst*(msp-1), but in general
!                                       will be much smaller (this parameter is similar to k13).
!                                 k19 - max number of compositions to be refined during adaptive
!                                       optimization.
!                                 k20 - max number of "z" coordinates saved for pseudocompounds
!                                       generated by adaptive refinement, max value = k21*m4, 
!                                       but usually much smaller (see k18).
!                                 k21 - max number of pseudocompounds for adaptive refinement.
!                                 k22 - 
!                                 k23 - max number of phases to be fractionated.
!                                 k24 - max number of discretization points for the standard simplex
!                                       in adaptive optimization.
      parameter (k0=25,k1=750000,k2=100000,k3=500,k4=23,k5=5)
      parameter (k7=k5+1,k8=k5+2) 
      parameter (k9=30,k10=240,k12=15,k13=mdim*k1,k14=18,k15=3,k16=30)
      parameter (k17=7,k18=29*k1)
      parameter (k19=2*k5+14,k21=1500000,k20=mdim*k21,k22=mdim*k19)
      parameter (k23=5,k24=1)  
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
      parameter (l2=5,l3=l2+2,l5=1000,l6=1000,l7=2048,l8=10)
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
      parameter (m0=8,m1=36,m2=8,m3=3,m4=15,m6=3,m7=12,m8=9,m9=10,m10=5)
      parameter (m11=5,m12=4,m13=8,m14=2,m15=9)
!                                 nsp - max number of species in fluid speciation routines 
      parameter (nsp=11)
!                                 nx - number of x-grid nodes in a contour data grid
!                                 ny - number of y-grid modes in a contour data grid
      parameter (nx=500,ny=500)
!
! NOTE: increasing parameter K5 requires changes to the following
! format statements:
!                     1000 in routine bulkck
!                     1000 in routine bulktst
!                     1150 in routine input1
!                     4000 in main in program build


!                                  the following are dependent parameters (also k7,k8)
      parameter (kd2=k8*31)

