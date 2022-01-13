!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_jpPotden
   ! Data type for the unsymmetrized potden variables. Needed for the port of
   ! juPhon (and applicable to SF magnetic fields)
   TYPE t_jpPotden
     INTEGER              :: vecdim ! 1 <-> scalar density/potential
                                    ! 3 <-> vectorial quantity
     INTEGER              :: dispatoms ! 1 <-> scalar quantity/gradient
                                        ! num_atoms <-> perturbed quantity
     COMPLEX, ALLOCATABLE :: pw(:,:,:,:), pw_w(:,:,:,:) ! iG, isp, idispat, idir
     COMPLEX, ALLOCATABLE :: mt(:,:,:,:,:,:) ! ir, ilm, iatom, isp, idispat, idir

   CONTAINS
     PROCEDURE :: init_jpPotden
     PROCEDURE :: resetjpPotden
     GENERIC   :: init=>init_jpPotden
     procedure :: jpdistribute
     procedure :: jpcollect
  END TYPE t_jpPotden

CONTAINS
  SUBROUTINE jpcollect(this,fmpi_comm)
    USE m_mpi_bc_tool
#ifdef CPP_MPI
    USE mpi
#endif
    IMPLICIT NONE
    CLASS(t_jpPotden), INTENT(INOUT) :: this
    INTEGER :: fmpi_comm
#ifdef CPP_MPI
    INTEGER :: ierr, irank
    COMPLEX, ALLOCATABLE :: ctmp(:)
    CALL MPI_COMM_RANK(fmpi_comm,irank,ierr)

    ALLOCATE(ctmp(size(this%pw)))
    CALL MPI_REDUCE(this%pw,ctmp,size(this%pw),MPI_DOUBLE_COMPLEX,MPI_SUM,0,fmpi_comm,ierr)
    IF (irank==0) this%pw=reshape(ctmp,shape(this%pw))
    DEALLOCATE(ctmp)

    ALLOCATE(ctmp(size(this%mt)))
    CALL MPI_REDUCE(this%mt,ctmp,size(this%mt),MPI_DOUBLE_COMPLEX,MPI_SUM,0,fmpi_comm,ierr)
    IF (irank==0) this%mt=reshape(ctmp,shape(this%mt))
    DEALLOCATE(ctmp)
#endif
END SUBROUTINE jpcollect

  SUBROUTINE jpdistribute(this,fmpi_comm)
    USE m_mpi_bc_tool
#ifdef CPP_MPI
    USE mpi
#endif
    IMPLICIT NONE
    CLASS(t_jpPotden), INTENT(INOUT) :: this
    INTEGER :: fmpi_comm
#ifdef CPP_MPI
    CALL mpi_bc(this%pw,0,fmpi_comm)
    IF (ALLOCATED(this%pw_w)) CALL mpi_bc(this%pw_w ,0,fmpi_comm)
    CALL mpi_bc(this%mt ,0,fmpi_comm)
#endif
END SUBROUTINE jpdistribute

  SUBROUTINE init_jpPotden(pd, vec_dim, disporder, nGq, jmtd, lmaxd, natoms, jspins, l_pw_w)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_jpPotden), INTENT(OUT) :: pd
    INTEGER, INTENT(IN)            :: vec_dim, disporder, nGq, jmtd, lmaxd, natoms, jspins
    LOGICAL,           INTENT(IN)  :: l_pw_w

    INTEGER:: err(4)

    err = 0
    pd%vecdim=vec_dim
    IF (disporder.EQ.0) THEN
      pd%dispatoms = 1
    else if (disporder.eq.1) then
      pd%dispatoms = natoms
    else
      CALL judft_error("Unreasonable density order!")
    END IF

    IF(ALLOCATED(pd%pw)) DEALLOCATE (pd%pw)
    IF(ALLOCATED(pd%pw_w)) DEALLOCATE (pd%pw_w)
    IF(ALLOCATED(pd%mt)) DEALLOCATE (pd%mt)

    ALLOCATE (pd%pw(nGq, jspins, pd%dispatoms, pd%vecdim),stat=err(1))
    IF (l_pw_w) THEN
        ALLOCATE (pd%pw_w(nGq, jspins, pd%dispatoms, pd%vecdim),stat=err(3))
    END IF
    ALLOCATE (pd%mt(jmtd, (lmaxd+1)**2, natoms, jspins, pd%dispatoms, pd%vecdim),stat=err(2))

    IF (ANY(err>0)) CALL judft_error("Not enough memory allocating potential or density")
    pd%pw=CMPLX(0.0,0.0)
    pd%mt=CMPLX(0.0,0.0)
  END SUBROUTINE init_jpPotden

  SUBROUTINE resetjpPotDen(pd)

    IMPLICIT NONE

    CLASS(t_jpPotden),INTENT(INOUT) :: pd

    pd%pw=CMPLX(0.0,0.0)
    pd%mt=0.0
    IF (ALLOCATED(pd%pw_w)) DEALLOCATE(pd%pw_w)
END SUBROUTINE resetjpPotDen

END MODULE m_types_jpPotden
