!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_fleurinput
  !this module 'uses' all types that extend t_fleurinput_base
  USE m_types_cell
  USE m_types_sym
  USE m_types_atoms
  USE m_types_input
  USE m_types_noco
  USE m_types_vacuum
  USE m_types_field
  USE m_types_sliceplot
  USE m_types_banddos
  USE m_types_mpinp
  USE m_types_hybinp
   
  USE m_types_coreSpecInput
  USE m_types_wannier
  USE m_types_xcpot
  USE m_types_forcetheo_data
  USE m_types_kpts
  USE m_types_enparaXML
  USE m_types_gfinp
  USE m_types_hub1inp
  USE m_types_juPhon
  IMPLICIT NONE

  TYPE t_fleurinput
    TYPE(t_cell)::cell
    TYPE(t_sym)::sym
    TYPE(t_atoms)::atoms
    TYPE(t_input)::input
    TYPE(t_noco)::noco
    TYPE(t_vacuum)::vacuum
    TYPE(t_field)::field
    TYPE(t_sliceplot)::sliceplot
    TYPE(t_banddos)::banddos
    TYPE(t_hybinp)::hybinp
    type(t_mpinp) :: mpinp
     
    TYPE(t_coreSpecInput)::coreSpecInput
    TYPE(t_forcetheo_data)::forcetheo_data
    TYPE(t_enparaXML)::enparaXML
    TYPE(t_kpts)::kpts
    type(t_gfinp)::gfinp
    type(t_hub1inp)::hub1inp
    type(t_juPhon)::juPhon
  end type t_fleurinput


END MODULE m_types_fleurinput
