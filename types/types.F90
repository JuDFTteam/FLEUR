!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!< This module collects all the type definitions. It can be used if no circular dependencies occur, i.e. in all
!! places not defining a type 
MODULE m_types
  USE m_types_rcmat
  USE m_types_xcpot
  USE m_types_lapw
  USE m_types_mpi
  USE m_types_tlmplm
  USE m_types_misc
  USE m_types_setup
  USE m_types_kpts
  USE m_types_usdus
  USE m_types_enpara
  USE m_types_potden
  USE m_types_cdnval
  USE m_types_field
  USE m_types_xcpot_inbuild
  USE m_types_regionCharges
  USE m_types_dos
  USE m_types_denCoeffsOffdiag
  USE m_types_force
  USE m_types_forcetheo
  USE m_types_hybrid
END MODULE m_types

