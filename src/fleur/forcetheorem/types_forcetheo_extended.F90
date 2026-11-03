!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_forcetheo_extended
  USE m_types_mae
  USE m_types_ssdisp
  USE m_types_dmi
  USE m_types_jij
  USE m_types_dmi_scf
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_forcetheo_mae, ssdisp_init, ssdisp_start, ssdisp_next_job, ssdisp_postprocess, ssdisp_dist, ssdisp_eval, &
     t_forcetheo_ssdisp, t_forcetheo_dmi, jij_init, jij_dist, jij_start, jij_next_job, jij_postprocess, jij_eval, &
     jij_q, map, fourier_transform, priv_analyse_data, t_forcetheo_jij, t_forcetheo_dmi_scf
END MODULE m_types_forcetheo_extended
