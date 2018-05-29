!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_diag
  USE m_juDFT
! the parameters are set to negative values to indicate that a particular solver is not compiled
    IMPLICIT NONE
    PRIVATE
#ifdef CPP_ELPA
  INTEGER,PARAMETER:: diag_elpa=1
#else
  INTEGER,PARAMETER:: diag_elpa=-1
#endif
#ifdef CPP_ELEMENTAL
  INTEGER,PARAMETER:: diag_elemental=2
#else
  INTEGER,PARAMETER:: diag_elemental=-2
#endif
#ifdef CPP_SCALAPACK
  INTEGER,PARAMETER:: diag_scalapack=3
#else
  INTEGER,PARAMETER:: diag_scalapack=-3
#endif
#ifdef CPP_MAGMA
  INTEGER,PARAMETER:: diag_magma=6
#else
  INTEGER,PARAMETER:: diag_magma=-6
#endif
  INTEGER,PARAMETER:: diag_lapack=4
  INTEGER,PARAMETER:: diag_chase=7
  INTEGER,PARAMETER:: diag_debugout=99
  PUBLIC eigen_diag,parallel_solver_available
CONTAINS

  LOGICAL FUNCTION parallel_solver_available()
    parallel_solver_available=any((/diag_elpa,diag_elemental,diag_scalapack/)>0)
  END FUNCTION parallel_solver_available

  SUBROUTINE eigen_diag(hmat,smat,ne,eig,ev)
    USE m_lapack_diag
    USE m_magma
    USE m_elpa
    USE m_scalapack
    USE m_elemental
    USE m_chase_diag
    USE m_types_mpimat
    IMPLICIT NONE
#ifdef CPP_MPI
    include 'mpif.h'
#endif
    CLASS(t_mat),INTENT(INOUT)             :: smat,hmat
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)   :: ev
    INTEGER,INTENT(INOUT)      :: ne
    REAL,INTENT(OUT)           :: eig(:)

    !Locals
    LOGICAL :: parallel
    SELECT TYPE(smat)
    CLASS IS (t_mpimat)
       parallel=.TRUE.
    CLASS default
       parallel=.FALSE.
    END SELECT
    
    CALL timestart("Diagonalization")
    !Select the solver
    SELECT CASE (priv_select_solver(parallel))
    CASE (diag_elpa)
       !CALL elpa(hmat,smat,ne,eig,ev)
    CASE (diag_elemental)
       !CALL ELEMENTAL(hmat,smat,ne,eig,ev)
    CASE (diag_scalapack)
       CALL scalapack(hmat,smat,ne,eig,ev)
    CASE (diag_magma)
       !CALL magma_diag(hmat,smat,ne,eig,ev)
    CASE (diag_lapack)
       CALL lapack_diag(hmat,smat,ne,eig,ev)
    CASE (diag_chase)
#ifdef CPP_CHASE
       CALL chase_diag(hmat,smat,ne,eig,ev)
#else
       CALL juDFT_error('ChASE eigensolver selected but not available', calledby = 'eigen_diag')
#endif
    CASE (diag_debugout)
       CALL priv_debug_out(smat,hmat)
    END SELECT
    CALL timestop("Diagonalization")
    !
  END SUBROUTINE eigen_diag


  SUBROUTINE priv_debug_out(smat,hmat)
    USE m_types
    use m_judft
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT) :: hmat,smat
    !small subroutine that does only wite the matrix to a file
    INTEGER:: i,ii,irank,ierr
    CHARACTER(len=20)::filename
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
#else
    irank=0
#endif
    WRITE(filename,"(a,i0)") "hmat",irank
    OPEN(999,file=TRIM(filename))
    WRITE(filename,"(a,i0)") "smat",irank
    OPEN(998,file=TRIM(filename))
    DO i=1,hmat%matsize2
       DO ii=1,hmat%matsize1
          IF (hmat%l_real) THEN
             WRITE(999,"(2i6,f15.6)") ii,i,hmat%data_r(ii,i)
             WRITE(998,"(2i6,f15.6)") ii,i,smat%data_r(ii,i)
          ELSE
             WRITE(999,"(2i6,2f15.6)") ii,i,hmat%data_c(ii,i)
             WRITE(998,"(2i6,2f15.6)") ii,i,smat%data_c(ii,i)
          ENDIF
       END DO
    ENDDO
    CLOSE(999)
    CLOSE(998)
#ifdef CPP_MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
    CALL judft_error("STOP in eigen_diag:debug_diag")
  END SUBROUTINE priv_debug_out

  FUNCTION priv_select_solver(parallel) RESULT(diag_solver)
    LOGICAL,INTENT(IN):: parallel
    INTEGER           :: diag_solver
    diag_solver=-99

    !Determine the default solver
    IF (parallel) THEN
#ifdef CPP_ELPA
       diag_solver=diag_elpa
#else
       diag_solver=diag_scalapack
#endif
    ELSE
#ifdef CPP_MAGMA
       diag_solver=diag_magma
#else
       diag_solver=diag_lapack
#endif
    ENDIF

    !check if a special solver was requested
    IF (juDFT_was_argument("-diag:elpa"))      diag_solver=diag_elpa
    IF (juDFT_was_argument("-diag:scalapack")) diag_solver=diag_scalapack
    IF (juDFT_was_argument("-diag:elemental")) diag_solver=diag_elemental
    IF (juDFT_was_argument("-diag:lapack"))    diag_solver=diag_lapack
    IF (juDFT_was_argument("-diag:magma"))     diag_solver=diag_magma
    IF (juDFT_was_argument("-diag:chase"))     diag_solver=diag_chase
    IF (juDFT_was_argument("-diag:debugout"))  diag_solver=diag_debugout
    
    !Check if solver is possible
    if (diag_solver<0) call priv_solver_error(diag_solver,parallel)
    IF (ANY((/diag_lapack,diag_magma/)==diag_solver).AND.parallel) CALL priv_solver_error(diag_solver,parallel)
    if (any((/diag_elpa,diag_elemental,diag_scalapack/)==diag_solver).and..not.parallel) call priv_solver_error(diag_solver,parallel)

  END FUNCTION priv_select_solver


  SUBROUTINE priv_solver_error(diag_solver,parallel)
    IMPLICIT NONE
    INTEGER,INTENT(IN):: diag_solver
    LOGICAL,INTENT(IN)::parallel
    SELECT CASE(diag_solver)
    CASE (diag_elpa)
       IF (parallel) THEN
          CALL juDFT_error("You did not compile with the ELPA solver and hence can not use it")
       ELSE
          CALL juDFT_error("The ELPA solver can not be used in serial")
       ENDIF
    CASE (diag_elemental)
       IF (parallel) THEN
          CALL juDFT_error("You did not compile with the ELEMENTAL solver and hence can not use it")
       ELSE
          CALL juDFT_error("The ELEMENTAL solver can not be used in serial")
       ENDIF
    CASE (diag_scalapack)
       IF (parallel) THEN
          CALL juDFT_error("You did not compile with the SCALAPACK solver and hence can not use it")
       ELSE
          CALL juDFT_error("The SCALAPACK solver can not be used in serial")
       ENDIF
    CASE (diag_lapack)
       IF (parallel) THEN
          CALL juDFT_error("The LAPACK solver can not be used in parallel")
       ENDIF
    CASE (diag_magma)
       CALL juDFT_error("You have not compiled with MAGMA support")
    CASE DEFAULT
       CALL judft_error("You have selected an unkown eigensolver")
    END SELECT
  END SUBROUTINE priv_solver_error


END MODULE m_eigen_diag
