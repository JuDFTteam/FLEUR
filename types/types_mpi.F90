!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mpi
  TYPE t_mpi
     INTEGER :: mpi_comm !< replaces MPI_COMM_WORLD
     INTEGER :: irank    !< rank of task in mpi_comm
     INTEGER :: isize    !< no of tasks in mpi_comm
     INTEGER :: n_start  !< no of first k-point to calculate on this PE
     INTEGER :: n_stride !< stride for k-loops
     INTEGER :: n_size   !< PE per kpoint, i.e. "isize" for eigenvalue parallelization
     INTEGER :: n_groups !< No of k-loops per PE
     INTEGER :: sub_comm !< Sub-Communicator for eigenvalue parallelization (all PE working on same k-point)
     INTEGER :: n_rank   !< rank in sub_comm
  END TYPE t_mpi
END MODULE m_types_mpi
