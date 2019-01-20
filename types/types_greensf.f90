
! This module defines a type for the Greens-functions used in the LDA+HIA formalism
! The Greens function is an on-Site Green' function which is stored in the matrix gmmpMat 
! and only contains Blocks with l = lprime in the MT-sphere

MODULE m_types_greensf

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensf 

         !Energy grid for Imaginary part
         INTEGER  :: ne       !number of energy grid points for imaginary part calculations
         REAL     :: e_top    !Cutoff energies
         REAL     :: e_bot
         REAL     :: sigma       !Smoothing parameter
     
         LOGICAL  :: l_tetra  !Determines wether to use the tetrahedron method for Brillouin-Zone integration

         !Energy contour parameters
         INTEGER  :: mode  !Determines the shape of the contour (more information in kkintgr.f90)
         INTEGER  :: nz    !number of points in the contour

         !array for energy contour
         COMPLEX, ALLOCATABLE  :: e(:)  !energy points
         COMPLEX, ALLOCATABLE  :: de(:) !weights for integration

         !Arrays for Green's function
         REAL, ALLOCATABLE :: im_gmmpMat(:,:,:,:,:)   !the imaginary part is stored in a different array because the number of energy points can differ
         COMPLEX, ALLOCATABLE :: gmmpMat(:,:,:,:,:)   
         REAL, ALLOCATABLE :: qalmmpMat(:,:,:,:,:,:) 

         CONTAINS
            PROCEDURE, PASS :: init => greensf_init
      END TYPE t_greensf

   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE greensf_init(thisGREENSF,input,atoms,kpts,dimension)

         USE m_types_setup
         USE m_types_kpts
         USE m_constants, only : lmaxU_const

         CLASS(t_greensf),       INTENT(INOUT)  :: thisGREENSF
         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_input),          INTENT(IN)     :: input
         TYPE(t_kpts),           INTENT(IN)     :: kpts
         TYPE(t_dimension),      INTENT(IN)     :: dimension

         INTEGER n

         !Parameters for calculation of the imaginary part
         thisGREENSF%ne       = input%ldahia_ne
         thisGREENSF%e_top    = input%ldahia_etop
         thisGREENSF%e_bot    = input%ldahia_ebot
         thisGREENSF%sigma    = input%ldahia_sigma

         thisGREENSF%l_tetra  = input%ldahia_tetra

         !Parameters for the energy contour in the complex plan
         !We use default values for now
         !thisGREENSF%mode     = input%ldahia_mode
         thisGREENSF%mode = 2
         n = 6

          IF(thisGREENSF%mode.EQ.1) THEN
            !thisGREENSF%nz = input%ldahia_nin
         ELSE IF(thisGREENSF%mode.EQ.2) THEN
            !n = input%ldahia_nin
            !ensure that we don't flood the memory accidentally
            IF(n.LT.2) n = 2
            !IF(n.GT.7) n = 7
            thisGREENSF%nz = 2**n
         END IF

         ALLOCATE (thisGREENSF%e(thisGREENSF%nz))
         ALLOCATE (thisGREENSF%de(thisGREENSF%nz))
         thisGREENSF%e(:) = CMPLX(0.0,0.0)
         thisGREENSF%de(:)= CMPLX(0.0,0.0)


         IF (thisGREENSF%l_tetra) THEN 
            ALLOCATE (thisGREENSF%qalmmpMat(dimension%neigd,kpts%nkpt,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
            thisGREENSF%qalmmpMat   = 0.0 
         ENDIF

         ALLOCATE (thisGREENSF%im_gmmpMat(thisGREENSF%ne,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
         ALLOCATE (thisGREENSF%gmmpMat(thisGREENSF%nz,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))

         thisGREENSF%im_gmmpMat     = 0.0
         thisGREENSF%gmmpMat     = CMPLX(0.0,0.0)

      END SUBROUTINE greensf_init

END MODULE m_types_greensf