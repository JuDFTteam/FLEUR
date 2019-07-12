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

         LOGICAL  :: l_calc   !Should the greens function be calculated in this iteration

         !Energy grid for Imaginary part
         INTEGER  :: ne       !number of energy grid points for imaginary part calculations
         REAL     :: e_top    !Cutoff energies
         REAL     :: e_bot
         REAL     :: del

         INTEGER, ALLOCATABLE :: kkintgr_cutoff(:,:,:)

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

      SUBROUTINE greensfCoeffs_init(thisGREENSFCOEFFS,input,lmax,atoms,noco,ef)

         USE m_juDFT
         USE m_types_setup

         CLASS(t_greensfCoeffs), INTENT(INOUT)  :: thisGREENSFCOEFFS
         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_input),          INTENT(IN)     :: input
         INTEGER,                INTENT(IN)     :: lmax
         REAL,                   INTENT(IN)     :: ef
         TYPE(t_noco),           INTENT(IN)     :: noco

         INTEGER i,j,l_dim,spin_dim

         !IF(thisGREENSFCOEFFS%l_calc) THEN 
         !
         !Set up general parameters for the Green's function (intersite and onsite)
         !
         !Parameters for calculation of the imaginary part
         thisGREENSFCOEFFS%ne       = input%gf_ne
         !take the energyParameterLimits from inp.xml if they are set, otherwise use default values
         IF(input%gf_ellow.NE.0.0.OR.input%gf_elup.NE.0.0) THEN
            thisGREENSFCOEFFS%e_top    = ef+input%gf_elup
            thisGREENSFCOEFFS%e_bot    = ef+input%gf_ellow
         ELSE
            thisGREENSFCOEFFS%e_top    = input%elup
            thisGREENSFCOEFFS%e_bot    = input%ellow
         ENDIF

         !set up energy grid for imaginary part
         thisGREENSFCOEFFS%del = (thisGREENSFCOEFFS%e_top-thisGREENSFCOEFFS%e_bot)/REAL(thisGREENSFCOEFFS%ne-1)

         spin_dim = MERGE(3,input%jspins,input%l_gfmperp)

         IF(atoms%n_gf.GT.0) THEN
            ALLOCATE(thisGREENSFCOEFFS%kkintgr_cutoff(atoms%n_gf,input%jspins,2))
            ALLOCATE (thisGREENSFCOEFFS%projdos(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim))
            thisGREENSFCOEFFS%projdos     = 0.0
            IF(.NOT.input%l_gfsphavg) THEN
               ALLOCATE (thisGREENSFCOEFFS%uu(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim))
               ALLOCATE (thisGREENSFCOEFFS%dd(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim))
               ALLOCATE (thisGREENSFCOEFFS%du(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim))
               ALLOCATE (thisGREENSFCOEFFS%ud(thisGREENSFCOEFFS%ne,MAX(1,atoms%n_gf),-lmax:lmax,-lmax:lmax,spin_dim))
               
               thisGREENSFCOEFFS%uu      = 0.0
               thisGREENSFCOEFFS%dd      = 0.0
               thisGREENSFCOEFFS%du      = 0.0  
               thisGREENSFCOEFFS%ud      = 0.0  
            ENDIF
         END IF

      END SUBROUTINE greensfCoeffs_init


END MODULE m_types_greensfCoeffs