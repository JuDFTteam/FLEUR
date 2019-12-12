!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_postprocessInput

CONTAINS

SUBROUTINE postprocessInput(mpi,input,field,sym,stars,atoms,vacuum,kpts,&
     oneD,cell,banddos,sliceplot,xcpot,forcetheo,forcetheo_data,&
     noco,sphhar,l_kpts)

  USE m_juDFT
  USE m_types
  USE m_constants
  USE m_ylm

  USE m_dwigner
  USE m_cdn_io
  USE m_prpxcfft


  use m_make_stars
  use m_make_sphhar


  USE m_convn
  USE m_efield
  USE m_od_kptsgen
  USE m_relaxio
  USE m_fleurinput_postprocess
  USE m_fleurinput_mpi_bc
  IMPLICIT NONE

  TYPE(t_mpi)      ,INTENT   (IN) :: mpi
  CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT):: forcetheo
  TYPE(t_forcetheo_data),INTENT(IN):: forcetheo_data
  TYPE(t_input),    INTENT(INOUT) :: input
  TYPE(t_sym),      INTENT(INOUT) :: sym
  TYPE(t_stars),    INTENT(INOUT) :: stars
  TYPE(t_atoms),    INTENT(INOUT) :: atoms
  TYPE(t_vacuum),   INTENT(INOUT) :: vacuum
  TYPE(t_kpts),     INTENT(INOUT) :: kpts
  TYPE(t_oneD),     INTENT(INOUT) :: oneD
  TYPE(t_cell),     INTENT(INOUT) :: cell
  TYPE(t_banddos),  INTENT(INOUT) :: banddos
  TYPE(t_sliceplot),INTENT(INOUT) :: sliceplot
  CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT) :: xcpot
  TYPE(t_noco),     INTENT(INOUT) :: noco

  TYPE(t_sphhar)   ,INTENT  (OUT) :: sphhar
  TYPE(t_field),    INTENT(INOUT) :: field
  LOGICAL,          INTENT   (IN) :: l_kpts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE postprocessInput

END MODULE m_postprocessInput
