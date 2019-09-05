!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_diag
  USE m_juDFT
  USE m_available_solvers
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eigen_diag

CONTAINS

  SUBROUTINE eigen_diag(solver,hmat,smat,ne,eig,ev,ikpt,jsp,iter)
    USE m_lapack_diag
    USE m_magma
    USE m_elpa
    USE m_elpa_onenode
    USE m_scalapack
    USE m_elemental
!    USE m_chase_diag
    USE m_types_mpimat
    USE m_types_gpumat
!    USE m_matrix_copy
    USE m_cusolver_diag
    USE m_judft_usage
    USE m_writeout
    IMPLICIT NONE
    INTEGER,                   INTENT(INOUT) :: solver
    CLASS(t_mat),              INTENT(INOUT) :: smat,hmat
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: ev         ! eigenvectors
    INTEGER,                   INTENT(INOUT) :: ne         ! number of eigenpairs to be found
    REAL,                      INTENT(OUT)   :: eig(:)     ! eigenvalues

    !Only for chase
    INTEGER,OPTIONAL,          INTENT(IN)    :: ikpt
    INTEGER,OPTIONAL,          INTENT(IN)    :: jsp
    INTEGER,OPTIONAL,          INTENT(IN)    :: iter

    !Locals  
    LOGICAL :: parallel

    SELECT TYPE(smat)
       CLASS IS (t_mpimat)
       parallel=.TRUE.
       CLASS default
       parallel=.FALSE.
    END SELECT

    solver=select_solver(solver,parallel)

    CALL timestart("Diagonalization")
    !Select the solver
    CALL add_usage_data("diag-solver", solver)
    SELECT CASE (solver)
    CASE (diag_elpa)
       CALL elpa_diag(hmat,smat,ne,eig,ev)
    CASE (diag_elpa_1node)
       CALL elpa_diag_onenode(hmat,smat,ne,eig,ev)
    CASE (diag_elemental)
       !CALL ELEMENTAL(hmat,smat,ne,eig,ev)
    CASE (diag_scalapack)
       CALL scalapack(hmat,smat,ne,eig,ev)
    CASE (diag_magma)
       CALL magma_diag(hmat,smat,ne,eig,ev)
    CASE (diag_cusolver)
       CALL cusolver_diag(hmat,smat,ne,eig,ev)
    CASE (diag_lapack)
       CALL lapack_diag(hmat,smat,ne,eig,ev)
    CASE (diag_chase)
       IF (.NOT.(PRESENT(ikpt).AND.PRESENT(jsp).AND.PRESENT(iter))) CALL judft_error("Optional arguments must be present for chase in eigen_diag")
!       CALL chase_diag(hmat,smat,ikpt,jsp,iter,ne,eig,ev)
    CASE (diag_debugout)
       CALL diag_writeout(smat,hmat)
    CASE default
       CALL judft_error("No solver available to diagonalize matrix")
    END SELECT
    CALL timestop("Diagonalization")

  END SUBROUTINE eigen_diag


END MODULE m_eigen_diag
