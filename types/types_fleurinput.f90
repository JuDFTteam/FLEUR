!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_fleurinput
   !*************************************************************
  !     This module contains the fleurinput-datatype.
  !     This type combines all setup-related data of FLEUR
  !     Usually this type should only be used for setting up a FLEUR calculation
   !*************************************************************
  USE m_types_cell
  USE m_types_sym
  USE m_types_banddos
  USE m_types_input
  USE m_types_sliceplot
  USE m_types_oneD
  USE m_types_hybrid
  USE m_types_noco
  USE m_types_stars
  USE m_types_atoms
  USE m_types_sphhar
  use m_types_dimension
  use m_types_coreSpecInput
  use m_types_wannier

  TODO:
  kpts
  xcpot
  forcetheo
  enpara
  
  implicit none
  private

  Type t_fleurinput
     type(t_cell)     ::cell
     type(t_sym)      ::sym
     type(t_atoms)    ::atoms
     type(t_input)    ::input
     
     type(t_banddos)  ::banddos
     type(t_sliceplot)::sliceplot
     type(t_oneD)     ::oneD
     type(t_hybrid)   ::hybrid
     type(t_noco)     ::noco
     type(t_stars)    ::stars
     type(t_sphhar)   ::sphhar
     type(t_dimension)::dimension
     type(t_wannier)  ::wannier
     type(t_coreSpecInput)::coreSpecInput
  end type t_fleurinput


end MODULE m_types_fleurinput
