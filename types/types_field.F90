!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_field
  !*************************************************************
  !     This module contains definitions for electric and magnetic field -types
  !*************************************************************
  PRIVATE
  REAL,TARGET::sigma=0.0
  TYPE t_efield
     REAL    :: zsigma  = 10.0  ! Distance to the charged plates
     REAL,POINTER    :: sigma   ! charge at the plates
     REAL    :: sig_b(2)=  0.0  ! Extra charge for the top/bottom plate
     COMPLEX :: vslope  =  0.0  ! Dirichlet bnd. cond.: Slope
     REAL,    ALLOCATABLE :: sigEF(:,:,:) ! (nx, ny, nvac)
     COMPLEX, ALLOCATABLE :: rhoEF(:,:)   ! (g_||, nvac)
     COMPLEX, ALLOCATABLE :: C1(:), C2(:) ! Coeff. for Dirichlet bnd.cond.
     LOGICAL :: l_segmented = .FALSE.
     LOGICAL :: plot_charge = .FALSE. ! Plot charge as inputted
     LOGICAL :: plot_rho    = .FALSE. ! Plot Fourier-transformed charge
     LOGICAL :: autocomp    = .TRUE.  ! Auto-compensate film charge
     LOGICAL :: dirichlet = .FALSE. ! Dirichlet vs. Neumann boundary cond.
     LOGICAL :: l_dirichlet_coeff = .FALSE. ! For MPI, true if C1/C2 set
  END TYPE t_efield

  TYPE t_field
     TYPE(t_efield)   :: efield
     LOGICAL          :: l_b_field=.false.
     REAL             :: b_field
     REAL,ALLOCATABLE :: b_field_mt(:)
   CONTAINS
     PROCEDURE :: init=>init_field
  END TYPE t_field

  PUBLIC t_field,t_efield
CONTAINS
  SUBROUTINE init_field(this,input)
    USE m_types_setup
    IMPLICIT NONE
    CLASS(t_field),INTENT(INOUT)::this
    TYPE(t_input),INTENT(INOUT) ::input
    input%sigma => sigma
    this%efield%sigma=>sigma
    PRINT *,"Sigma OK"
  END SUBROUTINE init_field
END MODULE m_types_field
