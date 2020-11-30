!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_xcpot_data
   !This module contains the xcpot-type used for the in-build xc-implementations
   IMPLICIT NONE

   TYPE t_xcpot_data
      !in the pbe case (exchpbe.F) lots of test are made
      !in addition some constants are set
      !to speed up this code precalculate things in init
      LOGICAL             :: is_rpbe = .false. !Rpbe
      LOGICAL             :: is_wc = .false.
      LOGICAL             :: is_hse = .false. !hse,lhse,vhse
      REAL                :: uk, um
      !many logicals to determine xcpot
      LOGICAL             :: is_pbes = .false.!is pbe-sol
      LOGICAL             :: is_pbe0 = .false.
      LOGICAL             :: is_bh = .false.
      LOGICAL             :: is_mjw = .false.
      REAL                :: exchange_weight
      INTEGER             :: krla !relativistic corrections
   contains 
      procedure :: mpi_bc => t_xcpot_data_mpi_bc 
   END TYPE t_xcpot_data
contains 

   subroutine t_xcpot_data_mpi_bc(data, rank, mpi_comm) 
      use m_mpi_bc_tool
      implicit NONE
      class(t_xcpot_data), intent(inout) :: data 
      integer, intent(in) :: rank, mpi_comm

      CALL mpi_bc(data%is_rpbe, rank, mpi_comm)
      CALL mpi_bc(data%is_wc, rank, mpi_comm)
      CALL mpi_bc(data%is_hse, rank, mpi_comm)
      CALL mpi_bc(data%uk, rank, mpi_comm)
      CALL mpi_bc(data%um, rank, mpi_comm)
      CALL mpi_bc(data%is_pbes, rank, mpi_comm)
      CALL mpi_bc(data%is_pbe0, rank, mpi_comm)
      CALL mpi_bc(data%is_bh, rank, mpi_comm)
      CALL mpi_bc(data%is_mjw, rank, mpi_comm)
      CALL mpi_bc(data%exchange_weight, rank, mpi_comm)
      CALL mpi_bc(data%krla, rank, mpi_comm)
   end subroutine t_xcpot_data_mpi_bc
END MODULE m_types_xcpot_data
