!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_potden
#ifdef CPP_MPI
   use mpi
#endif
CONTAINS
   SUBROUTINE mpi_bc_potden(fmpi,stars,sphhar,atoms,input,vacuum ,noco,potden,nococonv)

   USE m_types
   USE m_constants
   IMPLICIT NONE

   TYPE(t_mpi),INTENT(IN)        :: fmpi
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_stars),INTENT(IN)      :: stars
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_atoms),INTENT(IN)      :: atoms
   TYPE(t_noco),INTENT(IN)       :: noco
    
   TYPE(t_potden),INTENT(INOUT)  :: potden
   TYPE(t_nococonv),INTENT(INOUT),optional :: nococonv

   INTEGER :: n, ierr
   LOGICAL :: l_nocoAlloc, l_denMatAlloc, l_vaczAlloc, l_pw_wAlloc

#ifdef CPP_MPI

   CALL MPI_BCAST(potden%iter,1,MPI_INTEGER,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(potden%potdenType,1,MPI_INTEGER,0,fmpi%mpi_comm,ierr)

   l_nocoAlloc = .FALSE.
   l_denMatAlloc = .FALSE.
   l_vaczAlloc = .FALSE.
   l_pw_wAlloc = .FALSE.
   IF(fmpi%irank.EQ.0) THEN
      IF (ALLOCATED(potden%mmpMat)) l_denMatAlloc = .TRUE.
      IF (ALLOCATED(potden%vac)) l_vaczAlloc = .TRUE.
      IF (ALLOCATED(potden%pw_w)) l_pw_wAlloc = .TRUE.
   END IF
   CALL MPI_BCAST(l_nocoAlloc,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(l_denMatAlloc,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(l_vaczAlloc,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(l_pw_wAlloc,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
   IF((fmpi%irank.NE.0).AND.l_denMatAlloc) THEN
      IF(.NOT.ALLOCATED(potden%mmpMat)) THEN
         ALLOCATE(potDen%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                                MAX(1,atoms%n_u+atoms%n_hia),MERGE(3,input%jspins,any(noco%l_unrestrictMT).OR.noco%l_mperp)))
      END IF
   END IF

   n = stars%ng3 * SIZE(potden%pw,2)
   CALL MPI_BCAST(potden%pw,n,MPI_DOUBLE_COMPLEX,0,fmpi%mpi_comm,ierr)

   n = size(potden%mt)
   CALL MPI_BCAST(potden%mt,n,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)

   IF (l_pw_wAlloc) THEN
      n = stars%ng3 * SIZE(potden%pw_w,2)
      CALL MPI_BCAST(potden%pw_w,n,MPI_DOUBLE_COMPLEX,0,fmpi%mpi_comm,ierr)
   END IF

   IF (l_vaczAlloc) THEN
      n = vacuum%nmzd * stars%ng2 * 2 * SIZE(potden%vac,4)
      CALL MPI_BCAST(potden%vac,n,MPI_DOUBLE_COMPLEX,0,fmpi%mpi_comm,ierr)
   END IF

   IF (l_denMatAlloc) THEN
      n = SIZE(potden%mmpMat,1) * SIZE(potden%mmpMat,2) * SIZE(potden%mmpMat,3) * SIZE(potden%mmpMat,4)
      CALL MPI_BCAST(potden%mmpMat,n,MPI_DOUBLE_COMPLEX,0,fmpi%mpi_comm,ierr)
   END IF

   IF (present(nococonv)) THEN
     CALL MPI_BCAST(nococonv%alph, size(nococonv%alph), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
     CALL MPI_BCAST(nococonv%beta, size(nococonv%beta), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
     CALL MPI_BCAST(nococonv%alphprev, size(nococonv%alph), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
     CALL MPI_BCAST(nococonv%betaprev, size(nococonv%beta), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
   ENDIF
#endif

   END SUBROUTINE mpi_bc_potden
END MODULE m_mpi_bc_potden
