
! This module defines a type for the Greens-functions used in the LDA+HIA formalism
! The Greens function is an on-Site Green' function which is stored in the matrix gmmpMat 
! and only contains Blocks with l = lprime in the MT-sphere

MODULE m_types_greensf

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensf 

         INTEGER  :: ne       !number of energy grid points for imaginary part calculations

         !Cutoff parameters for energy integration
         REAL     :: e_top
         REAL     :: e_bot

         INTEGER  :: kkintgr_cut !cutoff for the kramers-kronig integration

         REAL     :: sigma       !Smoothing parameter

         LOGICAL  :: l_tetra  !Determines wether to use the tetrahedron method for Brillouin-Zone integration (not yet implemented)

         !Energy contour parameters

         INTEGER  :: mode 
         INTEGER  :: n_in

         INTEGER  :: nz

         COMPLEX, ALLOCATABLE  :: e(:)
         COMPLEX, ALLOCATABLE  :: de(:)



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

         !Cutoffs need to be calculated and sigma can probably taken from somewhere else
         thisGREENSF%ne       = input%ldahia_ne
         thisGREENSF%e_top    = input%ldahia_etop
         thisGREENSF%e_bot    = input%ldahia_ebot
         thisGREENSF%sigma    = input%ldahia_sigma
         thisGREENSF%l_tetra  = input%ldahia_tetra
         thisGREENSF%mode     = input%ldahia_mode
         thisGREENSF%n_in     = input%ldahia_nz


         ALLOCATE (thisGREENSF%gmmpMat(thisGREENSF%ne,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
         IF (thisGREENSF%l_tetra) THEN 
            ALLOCATE (thisGREENSF%qalmmpMat(dimension%neigd,kpts%nkpt,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
            thisGREENSF%qalmmpMat   = 0.0 
         ENDIF

         IF(thisGREENSF%mode.EQ.1) THEN
            ALLOCATE (thisGREENSF%e(thisGREENSF%n_in))
            ALLOCATE (thisGREENSF%de(thisGREENSF%n_in))

            thisGREENSF%e(:) = CMPLX(0.0,0.0)
            thisGREENSF%de(:)= CMPLX(0.0,0.0)
         ELSE IF(thisGREENSF%mode.EQ.2) THEN
            ALLOCATE (thisGREENSF%e(2**thisGREENSF%n_in))
            ALLOCATE (thisGREENSF%de(2**thisGREENSF%n_in))

            thisGREENSF%e(:) = CMPLX(0.0,0.0)
            thisGREENSF%de(:)= CMPLX(0.0,0.0)
         END IF


         thisGREENSF%gmmpMat     = CMPLX(0.0,0.0)
      END SUBROUTINE greensf_init

END MODULE m_types_greensf