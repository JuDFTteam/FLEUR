!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!< This module collects all the type definitions. It can be used if no circular dependencies occur, i.e. in all
!! places not defining a type
MODULE m_types
  USE m_types_cdnval
  USE m_types_denCoeffsOffdiag
  USE m_types_dos
  use m_types_eigvec
  USE m_types_enpara
  use m_types_fft
  USE m_types_fftGrid
  USE m_types_field
  USE m_types_force
  USE m_types_forcetheo
  USE m_types_greensf
  USE m_types_greensfCoeffs
  USE m_types_hub1data
  USE m_types_hybdat
  USE m_types_hybinp
  use m_types_hybmpi
  USE m_types_kpts
  USE m_types_lapw
  USE m_types_mat
  USE m_types_misc
  USE m_types_mpdata
  USE m_types_mpi
  USE m_types_mpinp
  USE m_types_nococonv
  USE m_types_potden
  USE m_types_regionCharges
  USE m_types_scalarGF
  USE m_types_setup
  USE m_types_tlmplm
  USE m_types_usdus
  USE m_types_xcpot
  USE m_types_xcpot_inbuild
END MODULE m_types
