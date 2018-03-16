!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_potden
CONTAINS
   SUBROUTINE mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,potden)

   USE m_types
   USE m_constants
   IMPLICIT NONE
   INCLUDE 'mpif.h'

   TYPE(t_mpi),INTENT(IN)        :: mpi
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_stars),INTENT(IN)      :: stars
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_atoms),INTENT(IN)      :: atoms
   TYPE(t_noco),INTENT(IN)       :: noco
   TYPE(t_oneD),INTENT(IN)       :: oneD
   TYPE(t_potden),INTENT(INOUT)  :: potden

   INTEGER :: n, ierr(3)
   LOGICAL :: l_nocoAlloc, l_denMatAlloc, l_vaczAlloc

   CALL MPI_BCAST(potden%iter,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)

   l_nocoAlloc = .FALSE.
   l_denMatAlloc = .FALSE.
   l_vaczAlloc = .FALSE.
   IF(mpi%irank.EQ.0) THEN
      IF (ALLOCATED(potden%cdom)) l_nocoAlloc = .TRUE.
      IF (ALLOCATED(potden%mmpMat)) l_denMatAlloc = .TRUE.
      IF (ALLOCATED(potden%vacz)) l_vaczAlloc = .TRUE.
   END IF
   CALL MPI_BCAST(l_nocoAlloc,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(l_denMatAlloc,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(l_vaczAlloc,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
   IF((mpi%irank.NE.0).AND.l_nocoAlloc) THEN
      IF (noco%l_noco) THEN
         IF(.NOT.ALLOCATED(potden%cdom)) ALLOCATE (potden%cdom(stars%ng3))
         IF(.NOT.ALLOCATED(potden%cdomvz)) ALLOCATE (potden%cdomvz(vacuum%nmzd,2))
      ELSE
         IF(.NOT.ALLOCATED(potden%cdom)) ALLOCATE (potden%cdom(1))
         IF(.NOT.ALLOCATED(potden%cdomvz)) ALLOCATE (potden%cdomvz(1,1))
      END IF
   END IF
   IF((mpi%irank.NE.0).AND.l_denMatAlloc) THEN
      IF(.NOT.ALLOCATED(potden%mmpMat)) THEN
         ALLOCATE(potDen%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_u),input%jspins))
      END IF
   END IF

   n = stars%ng3 * SIZE(potden%pw,2)
   CALL MPI_BCAST(potden%pw,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)

   n = atoms%jmtd * (sphhar%nlhd+1) * atoms%ntype * input%jspins
   CALL MPI_BCAST(potden%mt,n,MPI_DOUBLE,0,mpi%mpi_comm,ierr)

   IF (l_vaczAlloc) THEN
      n = vacuum%nmzd * 2 * SIZE(potden%vacz,3)
      CALL MPI_BCAST(potden%vacz,n,MPI_DOUBLE,0,mpi%mpi_comm,ierr)

      n = vacuum%nmzxyd * (stars%ng2-1) * 2 * SIZE(potden%vacxy,4)
      CALL MPI_BCAST(potden%vacxy,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
   END IF

   IF (l_nocoAlloc) THEN
      n = SIZE(potden%cdom,1)
      CALL MPI_BCAST(potden%cdom,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)

      n = SIZE(potden%cdomvz,1) * SIZE(potden%cdomvz,2)
      CALL MPI_BCAST(potden%cdomvz,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
   END IF

   IF (l_denMatAlloc) THEN
      n = SIZE(potden%mmpMat,1) * SIZE(potden%mmpMat,2) * SIZE(potden%mmpMat,3) * SIZE(potden%mmpMat,4)
      CALL MPI_BCAST(potden%mmpMat,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
   END IF

   END SUBROUTINE mpi_bc_potden
END MODULE m_mpi_bc_potden
