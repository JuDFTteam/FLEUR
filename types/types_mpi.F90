!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mpi
  TYPE t_mpi
     !k-point parallelism
     INTEGER :: mpi_comm !< replaces MPI_COMM_WORLD
     INTEGER :: irank    !< rank of task in mpi_comm
     INTEGER :: isize    !< no of tasks in mpi_comm
     INTEGER,ALLOCATABLE :: k_list(:)
     !Eigenvalue parallelism
     INTEGER :: sub_comm !< Sub-Communicator for eigenvalue parallelization (all PE working on same k-point)
     INTEGER :: n_rank   !< rank in sub_comm
     INTEGER :: n_size   !< PE per kpoint, i.e. "isize" for eigenvalue parallelization
     INTEGER,ALLOCATABLE :: ev_list(:)
  END TYPE t_mpi
END MODULE m_types_mpi
