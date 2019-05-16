MODULE m_types_greensfCoeffs

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_types_greensfCoeffs
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  Contains a type, which stores coefficients for the Green's function calculated 
   !>  in the k-point loop in cdnval
   !>  Contains Arrays for the following cases:
   !>       -onsite
   !>           -spherically averaged/radial dependence (r=r')
   !>           -non-magnetic/collinear/noco
   !>       -intersite
   !>  Furthermore this module contains the information about the energy grid where
   !>  the imaginary part is calculated
   ! REVISION HISTORY:
   ! February 2019 - Initial Version
   !------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensfCoeffs

         LOGICAL  :: l_onsite !Are the arrays for onsite calculations allocated
         LOGICAL  :: l_intersite !Are the arrays for intersite calculations allocated

         !Energy grid for Imaginary part
         INTEGER  :: ne       !number of energy grid points for imaginary part calculations
         REAL     :: e_top    !Cutoff energies
         REAL     :: e_bot
         REAL     :: del
         REAL     :: sigma       !Smoothing parameter(not used at the moment)

         INTEGER, ALLOCATABLE :: kkintgr_cutoff(:,:)

         !Array declarations 
         !If we look at the Green's function that only depends on Energy and not on spatial arguments
         !the imaginary part is equal to the proected density of states
         REAL, ALLOCATABLE :: projdos(:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         REAL, ALLOCATABLE :: uu(:,:,:,:,:)
         REAL, ALLOCATABLE :: dd(:,:,:,:,:)
         REAL, ALLOCATABLE :: du(:,:,:,:,:)
         REAL, ALLOCATABLE :: ud(:,:,:,:,:)

         CONTAINS
            PROCEDURE, PASS :: init => greensfCoeffs_init
      END TYPE t_greensfCoeffs

   PUBLIC t_greensfCoeffs

   CONTAINS

      SUBROUTINE greensfCoeffs_init(thisGREENSFCOEFFS,input,lmax,atoms,kpts,noco,l_onsite,l_intersite)

         USE m_juDFT
         USE m_types_setup
         USE m_types_kpts
         USE m_constants, only : lmaxU_const

         CLASS(t_greensfCoeffs), INTENT(INOUT)  :: thisGREENSFCOEFFS
         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_input),          INTENT(IN)     :: input
         INTEGER,                INTENT(IN)     :: lmax
         TYPE(t_kpts), OPTIONAL, INTENT(IN)     :: kpts
         TYPE(t_noco), OPTIONAL, INTENT(IN)     :: noco
         LOGICAL,                INTENT(IN)     :: l_onsite
         LOGICAL,                INTENT(IN)     :: l_intersite

         INTEGER i,j,l_dim
         REAL    tol
         LOGICAL l_new
         INTEGER n_spins

         thisGREENSFCOEFFS%l_onsite = l_onsite
         thisGREENSFCOEFFS%l_intersite = l_intersite

         IF(.NOT.l_onsite.AND.noco%l_mperp) CALL juDFT_error("NOCO + intersite gf not implented",calledby="greensf_init")
         IF(l_onsite.AND.noco%l_mperp) CALL juDFT_error("NOCO + onsite gf not implented",calledby="greensf_init")

         !
         !Set up general parameters for the Green's function (intersite and onsite)
         !
         tol = 1e-14
         !Parameters for calculation of the imaginary part
         thisGREENSFCOEFFS%ne       = input%onsite_ne
         !take the energyParameterLimits from inp.xml if they are set, otherwise use default values
         IF(ABS(input%ellow).LT.tol.AND.ABS(input%elup).LT.tol) THEN
            thisGREENSFCOEFFS%e_top    = 2.0
            thisGREENSFCOEFFS%e_bot    = -2.0
         ELSE
            thisGREENSFCOEFFS%e_top    = input%elup
            thisGREENSFCOEFFS%e_bot    = input%ellow
         ENDIF
         thisGREENSFCOEFFS%e_top    = 2.0
         thisGREENSFCOEFFS%e_bot    = -2.0
         thisGREENSFCOEFFS%sigma    = input%onsite_sigma

         !set up energy grid for imaginary part
         thisGREENSFCOEFFS%del      = (thisGREENSFCOEFFS%e_top-thisGREENSFCOEFFS%e_bot)/REAL(thisGREENSFCOEFFS%ne-1)

         IF(thisGREENSFCOEFFS%l_onsite) THEN
            !
            !In the case of an onsite gf we look at the case l=l' and r=r' on one site
            !

            IF(atoms%n_gf.GT.0) THEN !Are there Green's functions to be calculated?
               !Do we need the off-diagonal elements
               ALLOCATE(thisGREENSFCOEFFS%kkintgr_cutoff(atoms%n_gf,2))
               n_spins = MERGE(4,input%jspins,noco%l_mperp)
               ALLOCATE (thisGREENSFCOEFFS%projdos(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_spins))
               thisGREENSFCOEFFS%projdos     = 0.0
               IF(.NOT.input%onsite_sphavg) THEN
                  ALLOCATE (thisGREENSFCOEFFS%uu(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_spins))
                  ALLOCATE (thisGREENSFCOEFFS%dd(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_spins))
                  ALLOCATE (thisGREENSFCOEFFS%du(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_spins))
                  ALLOCATE (thisGREENSFCOEFFS%ud(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_spins))
                  
                  thisGREENSFCOEFFS%uu      = 0.0
                  thisGREENSFCOEFFS%dd      = 0.0
                  thisGREENSFCOEFFS%du      = 0.0  
                  thisGREENSFCOEFFS%ud      = 0.0  
               ENDIF
            ENDIF
         END IF

         IF(thisGREENSFCOEFFS%l_intersite) THEN
            !
            !intersite case; Here we look at l/=l' r/=r' and multiple sites
            !

            l_dim = lmax**2 + lmax
            !How to incorporate smoothly into the arrays

         ENDIF 

      END SUBROUTINE greensfCoeffs_init


END MODULE m_types_greensfCoeffs