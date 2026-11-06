!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_types_desym_map
   use m_mpi_bc_tool
   implicit none
   private
   public :: t_desym_map

   ! Relation between a symmetric setup and its copy without symmetry. Atom indices are identical in both setups.
   type t_desym_map
      integer, allocatable :: parent_type(:)  ! symmetric type of each nosym type (= atom)
      integer, allocatable :: type_op(:)      ! operation mapping the representative atom onto each atom (sym%ngopr)
      integer, allocatable :: ldau_parent(:)  ! symmetric LDA+U entry of each nosym LDA+U entry
      integer, allocatable :: kpt_parent(:)   ! IBZ k index of each full-BZ k (kpts%bkp)
      integer, allocatable :: kpt_op(:)       ! operation generating it, > sym%nop means time reversal (kpts%bksym)
   contains
      procedure :: mpi_bc => mpi_bc_desym_map
   end type t_desym_map

contains

   subroutine mpi_bc_desym_map(this, mpi_comm, irank)
      !! Broadcasts all maps from rank irank (default 0).
      class(t_desym_map), intent(inout) :: this
      integer,            intent(in)    :: mpi_comm
      integer, optional,  intent(in)    :: irank

      integer :: rank

      rank = 0
      if (present(irank)) rank = irank

      call mpi_bc(this%parent_type, rank, mpi_comm)
      call mpi_bc(this%type_op, rank, mpi_comm)
      call mpi_bc(this%ldau_parent, rank, mpi_comm)
      call mpi_bc(this%kpt_parent, rank, mpi_comm)
      call mpi_bc(this%kpt_op, rank, mpi_comm)

   end subroutine mpi_bc_desym_map

end module m_types_desym_map
