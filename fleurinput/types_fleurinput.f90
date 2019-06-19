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
  USE m_types_hybrid
  USE m_types_oneD
  USE m_types_coreSpecInput
  USE m_types_wannier
  USE m_types_xcpot
  USE m_types_forcetheo_data
  USE m_types_kpts
  USE m_types_enparaXML
  IMPLICIT NONE

CONTAINS
  !Subroutine does nothing, only here for copy-paste code...
  SUBROUTINE dummy(cell,sym,atoms,input,noco,vacuum,field,&
       sliceplot,banddos,hybrid,oneD,coreSpecInput,wann,&
       xcpot,forcetheo_data,kpts,enparaXML)
    TYPE(t_cell),INTENT(IN)::cell
    TYPE(t_sym),INTENT(IN)::sym
    TYPE(t_atoms),INTENT(IN)::atoms
    TYPE(t_input),INTENT(IN)::input
    TYPE(t_noco),INTENT(IN)::noco
    TYPE(t_vacuum),INTENT(IN)::vacuum
    TYPE(t_field),INTENT(IN)::field
    TYPE(t_sliceplot),INTENT(IN)::sliceplot
    TYPE(t_banddos),INTENT(IN)::banddos
    TYPE(t_hybrid),INTENT(IN)::hybrid
    TYPE(t_oneD),INTENT(IN)::oneD
    TYPE(t_coreSpecInput),INTENT(IN)::coreSpecInput
    TYPE(t_wann),INTENT(IN)::wann
    CLASS(t_xcpot),INTENT(IN)::xcpot
    TYPE(t_forcetheo_data),INTENT(IN)::forcetheo_data
    TYPE(t_enparaXML),INTENT(IN)::enparaXML
    TYPE(t_kpts),INTENT(IN)::kpts
  END SUBROUTINE dummy
END MODULE m_types_fleurinput

