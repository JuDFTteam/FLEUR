MODULE m_types_greensfCoeffs

   !This type contains all the information used in the calculation of the green's function from the DFT eigenstates
   !These arrays were split off from t_greensf to avoid dragging these temporary arrays through the whole program

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
     
         LOGICAL  :: l_tetra  !Determines wether to use the tetrahedron method for Brillouin-Zone integration

         !we store the atom types and l's for which to calculate the onsite gf to make it easier to reuse in other circumstances
         INTEGER  :: n_gf
         INTEGER, ALLOCATABLE :: atomType(:)
         INTEGER, ALLOCATABLE :: l_gf(:)

         !Array declarations
         REAL, ALLOCATABLE :: im_g(:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         REAL, ALLOCATABLE :: uu(:,:,:,:,:)
         REAL, ALLOCATABLE :: dd(:,:,:,:,:)
         REAL, ALLOCATABLE :: du(:,:,:,:,:)
         REAL, ALLOCATABLE :: ud(:,:,:,:,:)

         !noco case
         REAL, ALLOCATABLE :: uu21(:,:,:,:)
         REAL, ALLOCATABLE :: dd21(:,:,:,:)
         REAL, ALLOCATABLE :: du21(:,:,:,:)
         COMPLEX, ALLOCATABLE :: im_g21(:,:,:,:)

         !intersite case
         REAL, ALLOCATABLE :: uu_int(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: dd_int(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: du_int(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: ud_int(:,:,:,:,:,:)

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
            thisGREENSFCOEFFS%e_top    = 1.0
            thisGREENSFCOEFFS%e_bot    = -1.0
         ELSE
            thisGREENSFCOEFFS%e_top    = input%elup
            thisGREENSFCOEFFS%e_bot    = input%ellow
         ENDIF

         thisGREENSFCOEFFS%sigma    = input%onsite_sigma
         thisGREENSFCOEFFS%l_tetra  = input%onsite_tetra

         !set up energy grid for imaginary part
         thisGREENSFCOEFFS%del      = (thisGREENSFCOEFFS%e_top-thisGREENSFCOEFFS%e_bot)/REAL(thisGREENSFCOEFFS%ne-1)

         !Determine for which types and l's to calculate the onsite gf
         ALLOCATE(thisGREENSFCOEFFS%atomType(MAX(1,atoms%n_hia+atoms%n_j0)))
         ALLOCATE(thisGREENSFCOEFFS%l_gf(MAX(1,atoms%n_hia+atoms%n_j0)))
         thisGREENSFCOEFFS%atomType(:) = 0
         thisGREENSFCOEFFS%l_gf(:) = 0

         IF(thisGREENSFCOEFFS%l_onsite) THEN
            !
            !In the case of an onsite gf we look at the case l=l' and r=r' on one site
            !
            thisGREENSFCOEFFS%n_gf = 0
            !DFT+HIA:
            DO i = 1, atoms%n_hia

               thisGREENSFCOEFFS%n_gf = thisGREENSFCOEFFS%n_gf + 1
               thisGREENSFCOEFFS%atomType(thisGREENSFCOEFFS%n_gf) =  atoms%lda_hia(i)%atomType
               thisGREENSFCOEFFS%l_gf(thisGREENSFCOEFFS%n_gf)     =  atoms%lda_hia(i)%l

            ENDDO

            !Effective exchange interaction:
            DO i = 1, atoms%n_j0
               !Avoid double calculations:
               l_new = .true.
               DO j = 1, thisGREENSFCOEFFS%n_gf
                  IF(thisGREENSFCOEFFS%atomType(j).EQ.atoms%j0(i)%atomType.AND.thisGREENSFCOEFFS%l_gf(j).EQ.atoms%j0(i)%l) THEN
                     l_new = .false.
                     EXIT
                  ENDIF
               ENDDO

               IF(l_new) THEN
                  thisGREENSFCOEFFS%n_gf = thisGREENSFCOEFFS%n_gf + 1
                  thisGREENSFCOEFFS%atomType(thisGREENSFCOEFFS%n_gf) =  atoms%j0(i)%atomType
                  thisGREENSFCOEFFS%l_gf(thisGREENSFCOEFFS%n_gf)     =  atoms%j0(i)%l
               ENDIF 
            ENDDO

            IF(thisGREENSFCOEFFS%n_gf.GT.0) THEN !Are there Green's functions to be calculated?

               IF(input%onsite_sphavg) THEN
                  ALLOCATE (thisGREENSFCOEFFS%im_g(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
                  thisGREENSFCOEFFS%im_g     = 0.0
               ELSE
                  ALLOCATE (thisGREENSFCOEFFS%uu(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
                  ALLOCATE (thisGREENSFCOEFFS%dd(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
                  ALLOCATE (thisGREENSFCOEFFS%du(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
                  
                  thisGREENSFCOEFFS%uu      = 0.0
                  thisGREENSFCOEFFS%dd      = 0.0
                  thisGREENSFCOEFFS%du      = 0.0  
               ENDIF

               !Allocate arrays for non-colinear part
               IF(noco%l_mperp) THEN
                  IF(.NOT.input%onsite_sphavg) THEN 
                     ALLOCATE (thisGREENSFCOEFFS%uu21(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
                     ALLOCATE (thisGREENSFCOEFFS%dd21(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
                     ALLOCATE (thisGREENSFCOEFFS%du21(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
                    
                     thisGREENSFCOEFFS%uu21    = 0.0
                     thisGREENSFCOEFFS%dd21    = 0.0
                     thisGREENSFCOEFFS%du21    = 0.0
                  ENDIF

                  ALLOCATE (thisGREENSFCOEFFS%im_g21(thisGREENSFCOEFFS%ne,MAX(1,thisGREENSFCOEFFS%n_gf),&
                                                            -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
                  thisGREENSFCOEFFS%im_g21  = 0.0
               ENDIF
            ENDIF
         END IF

         IF(thisGREENSFCOEFFS%l_intersite) THEN
            !
            !intersite case; Here we look at l/=l' r/=r' and multiple sites
            !

            l_dim = lmax**2 + lmax
            ALLOCATE (thisGREENSFCOEFFS%uu_int(thisGREENSFCOEFFS%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))
            ALLOCATE (thisGREENSFCOEFFS%dd_int(thisGREENSFCOEFFS%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))
            ALLOCATE (thisGREENSFCOEFFS%du_int(thisGREENSFCOEFFS%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))
            ALLOCATE (thisGREENSFCOEFFS%ud_int(thisGREENSFCOEFFS%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))

            thisGREENSFCOEFFS%uu_int      = 0.0
            thisGREENSFCOEFFS%dd_int      = 0.0
            thisGREENSFCOEFFS%du_int      = 0.0
            thisGREENSFCOEFFS%ud_int      = 0.0

         ENDIF 

      END SUBROUTINE greensfCoeffs_init


END MODULE m_types_greensfCoeffs