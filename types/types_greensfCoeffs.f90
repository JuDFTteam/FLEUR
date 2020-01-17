!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

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

   USE m_juDFT
   USE m_types_setup
   USE m_constants

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensfCoeffs

         INTEGER, ALLOCATABLE :: kkintgr_cutoff(:,:,:)

         !Array declarations
         !If we look at the Green's function that only depends on Energy and not on spatial arguments
         !the imaginary part is equal to the proected density of states
         COMPLEX, ALLOCATABLE :: projdos(:,:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: du(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:,:)

         CONTAINS
            PROCEDURE, PASS :: init    =>  greensfCoeffs_init
      END TYPE t_greensfCoeffs

   PUBLIC t_greensfCoeffs

   CONTAINS

      SUBROUTINE greensfCoeffs_init(thisGREENSFCOEFFS,gfinp,input,atoms,noco)

         CLASS(t_greensfCoeffs), INTENT(INOUT)  :: thisGREENSFCOEFFS
         TYPE(t_gfinp),          INTENT(IN)     :: gfinp
         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_input),          INTENT(IN)     :: input
         TYPE(t_noco),           INTENT(IN)     :: noco

         INTEGER lmax,spin_dim

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const

         IF(gfinp%n.GT.0) THEN
            ALLOCATE(thisGREENSFCOEFFS%kkintgr_cutoff(gfinp%n,input%jspins,2),source=0)
            ALLOCATE (thisGREENSFCOEFFS%projdos(gfinp%ne,-lmax:lmax,-lmax:lmax,0:MAXVAL(atoms%neq),MAX(1,gfinp%n),spin_dim),source=cmplx_0)
            IF(.NOT.input%l_gfsphavg) THEN
               ALLOCATE (thisGREENSFCOEFFS%uu(gfinp%ne,-lmax:lmax,-lmax:lmax,0:MAXVAL(atoms%neq),MAX(1,gfinp%n),spin_dim),source=cmplx_0)
               ALLOCATE (thisGREENSFCOEFFS%dd(gfinp%ne,-lmax:lmax,-lmax:lmax,0:MAXVAL(atoms%neq),MAX(1,gfinp%n),spin_dim),source=cmplx_0)
               ALLOCATE (thisGREENSFCOEFFS%du(gfinp%ne,-lmax:lmax,-lmax:lmax,0:MAXVAL(atoms%neq),MAX(1,gfinp%n),spin_dim),source=cmplx_0)
               ALLOCATE (thisGREENSFCOEFFS%ud(gfinp%ne,-lmax:lmax,-lmax:lmax,0:MAXVAL(atoms%neq),MAX(1,gfinp%n),spin_dim),source=cmplx_0)
            ENDIF
         END IF

      END SUBROUTINE greensfCoeffs_init

END MODULE m_types_greensfCoeffs