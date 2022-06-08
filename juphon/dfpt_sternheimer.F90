!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_sternheimer
   USE m_types
   USE m_make_stars
   USE m_vgen
   USE m_eigen
   USE m_dfpt_cdngen
   USE m_mix

IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt_sternheimer(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, &
                               rho, vTot, grRho, grVtot, q_list, iQ, iType, iDir, &
                               rho1, vTot1, rho1Im, vTot1Im)
      TYPE(t_fleurinput), INTENT(IN)    :: fi
      TYPE(t_sphhar),     INTENT(IN)    :: sphhar
      TYPE(t_stars),      INTENT(IN)    :: stars
      TYPE(t_nococonv),   INTENT(IN)    :: nococonv
      TYPE(t_kpts),       INTENT(IN)    :: qpts
      TYPE(t_mpi),        INTENT(IN)    :: fmpi
      TYPE(t_results),    INTENT(INOUT) :: results
      TYPE(t_enpara),     INTENT(IN)    :: enpara
      TYPE(t_potden),     INTENT(IN)    :: rho, vTot, grRho, grVtot
      TYPE(t_potden),     INTENT(INOUT) :: rho1, vTot1, rho1Im, vTot1Im

      INTEGER, INTENT(IN) :: q_list(:)
      INTEGER, INTENT(IN) :: iQ, iType, iDir

      TYPE(t_stars) :: starsq

      CALL make_stars(starsq, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qpts%bk(:,iQ))
      CALL rho1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
      CALL rho1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)

   END SUBROUTINE dfpt_sternheimer
END MODULE m_dfpt_sternheimer
