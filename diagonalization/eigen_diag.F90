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
#ifdef CPP_CHASE
  INTEGER,PARAMETER:: diag_chase=7
#else
  INTEGER,PARAMETER:: diag_chase=-7
#endif
#ifdef CPP_CUSOLVER
  INTEGER,PARAMETER:: diag_cusolver=8
#else
  INTEGER,PARAMETER:: diag_cusolver=-8
#endif
  
  INTEGER,PARAMETER:: diag_lapack=4
  INTEGER,PARAMETER:: diag_elpa_1node=14
 
  INTEGER,PARAMETER:: diag_debugout=99
  PUBLIC eigen_diag,parallel_solver_available
CONTAINS

  LOGICAL FUNCTION parallel_solver_available()
    parallel_solver_available=any((/diag_elpa,diag_elemental,diag_scalapack/)>0)
  END FUNCTION parallel_solver_available

  SUBROUTINE eigen_diag(mpi,hmat,smat,ikpt,jsp,iter,ne,eig,ev)
    USE m_lapack_diag
    USE m_magma
    USE m_elpa
    USE m_elpa_onenode
    USE m_scalapack
    USE m_elemental
    USE m_chase_diag
    USE m_types_mpimat
    USE m_types_gpumat
    USE m_matrix_copy
    USE m_cusolver_diag
    IMPLICIT NONE
#ifdef CPP_MPI
    include 'mpif.h'
#endif
    TYPE(t_mpi),               INTENT(IN)    :: mpi
    CLASS(t_mat),              INTENT(INOUT) :: smat,hmat
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: ev
    INTEGER,                   INTENT(IN)    :: ikpt
    INTEGER,                   INTENT(IN)    :: jsp
    INTEGER,                   INTENT(IN)    :: iter
    INTEGER,                   INTENT(INOUT) :: ne
    REAL,                      INTENT(OUT)   :: eig(:)

    !Locals  
    LOGICAL :: parallel
    !For check-mode
    TYPE(t_mat) :: s_store,h_store
    
         
    SELECT TYPE(smat)
    CLASS IS (t_mpimat)
       parallel=.TRUE.
    CLASS default
       parallel=.FALSE.
    END SELECT

    !Create a copy of the matrix if in test mode
    IF (TRIM(judft_string_for_argument("-diag"))=="test") THEN
       SELECT TYPE(hmat)
        CLASS IS (t_mpimat)
          CALL s_store%init(hmat%l_real,hmat%global_size1,hmat%global_size2)
          CALL h_store%init(hmat%l_real,hmat%global_size1,hmat%global_size2)
       CLASS default
          CALL s_store%init(hmat%l_real,hmat%matsize1,hmat%matsize2)
          CALL h_store%init(hmat%l_real,hmat%matsize1,hmat%matsize2)
       END SELECT
       CALL matrix_copy(smat,s_store)
       CALL matrix_copy(hmat,h_store)
    END IF
       
    CALL timestart("Diagonalization")
    !Select the solver
    SELECT CASE (priv_select_solver(parallel))
    CASE (diag_elpa)
       CALL elpa_diag(hmat,smat,ne,eig,ev)
    CASE (diag_elpa_1node)
       CALL elpa_diag_onenode(hmat,smat,ne,eig,ev)
    CASE (diag_elemental)
       !CALL ELEMENTAL(hmat,smat,ne,eig,ev)
    CASE (diag_scalapack)
       CALL scalapack(hmat,smat,ne,eig,ev)
    CASE (diag_magma)
       !CALL magma_diag(hmat,smat,ne,eig,ev)
    CASE (diag_cusolver)
       CALL cusolver_diag(hmat,smat,ne,eig,ev)
    CASE (diag_lapack)
       CALL lapack_diag(hmat,smat,ne,eig,ev)
    CASE (diag_chase)
       CALL chase_diag(hmat,smat,ikpt,jsp,iter,ne,eig,ev)
    CASE (diag_debugout)
       CALL priv_debug_out(smat,hmat)
    END SELECT

    !Create test the solution
    IF (TRIM(judft_string_for_argument("-diag"))=="test") THEN
       CALL priv_test_solution(mpi,ne,s_store,h_store,eig,ev)
       CALL judft_error("Diagonalization tested")
    END IF
    CALL timestop("Diagonalization")
    !
  END SUBROUTINE eigen_diag


  SUBROUTINE priv_test_solution(mpi,ne,s_store,h_store,eig1,ev)
    USE m_types
    USE m_lapack_diag
    USE m_matrix_copy
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(in):: mpi
    INTEGER,INTENT(INOUT) :: ne
    TYPE(t_mat)           :: s_store,h_store
    REAL,INTENT(in)       :: eig1(:)
    CLASS(t_mat)          :: ev

    
    REAL,ALLOCATABLE::eig2(:)
    TYPE(t_mat)     ::ev1
    CLASS(t_mat),ALLOCATABLE    ::ev2
    INTEGER         :: i,irank
    

    ALLOCATE(eig2(ne))
    CALL ev1%init(s_store%l_real,s_store%matsize1,s_store%matsize2)
    CALL matrix_copy(ev,ev1)
    
    IF (mpi%irank==0) THEN
       CALL lapack_diag(h_store,s_store,ne,eig2,ev2)

       OPEN(99,file="diag.compare")
       WRITE(99,*) "Eigenvalues"
       DO i=1,ne
          WRITE(99,*) i,eig1(i),eig2(i)
       ENDDO
       WRITE(99,*) "Eigenvectors"
       DO i=1,ne
          IF (ev1%l_real) THEN
             WRITE(99,"(i0,20(1x,f10.5))") i,ev1%data_r(1:10,i)
             WRITE(99,"(i0,20(1x,f10.5))") i,ev2%data_r(1:10,i)
          ELSE
             WRITE(99,"(i0,20(1x,f10.5))") i,ev1%data_c(1:10,i)
             WRITE(99,"(i0,20(1x,f10.5))") i,ev2%data_c(1:10,i)
          END IF
       ENDDO
       CLOSE(99)
    END IF
  END SUBROUTINE priv_test_solution
  
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
       diag_solver=diag_lapack
    ENDIF

    !check if a special solver was requested
    IF (TRIM(juDFT_string_for_argument("-diag"))=="elpa")       diag_solver=diag_elpa
    IF (TRIM(juDFT_string_for_argument("-diag"))=="elpa_1node") diag_solver=diag_elpa_1node
    IF (trim(juDFT_string_for_argument("-diag"))=="scalapack")  diag_solver=diag_scalapack
    IF (trim(juDFT_string_for_argument("-diag"))=="elemental")  diag_solver=diag_elemental
    IF (trim(juDFT_string_for_argument("-diag"))=="lapack")     diag_solver=diag_lapack
    IF (trim(juDFT_string_for_argument("-diag"))=="magma")      diag_solver=diag_magma
    IF (trim(juDFT_string_for_argument("-diag"))=="chase")      diag_solver=diag_chase
    IF (trim(juDFT_string_for_argument("-diag"))=="cusolver")   diag_solver=diag_cusolver
    IF (trim(juDFT_string_for_argument("-diag"))=="debugout")   diag_solver=diag_debugout
    
    !Check if solver is possible
    IF (diag_solver<0)  CALL juDFT_error("You selected a solver for the eigenvalue problem that is not available",hint="You most probably did not provide the appropriate libraries for compilation/linking")
    IF (ANY((/diag_lapack,diag_magma,diag_cusolver/)==diag_solver).AND.parallel) CALL judft_error("You selected an eigensolver that does not support distributed memory parallism",hint="Try scalapack,elpa or another supported solver for parallel matrices")
    IF (ANY((/diag_elpa,diag_elemental,diag_scalapack/)==diag_solver).AND..NOT.parallel) CALL judft_error("You selected an eigensolver for matrices that are memory distributed",hint="Try lapack, cusolver or another supported solver for non-distributed matrices")


  END FUNCTION priv_select_solver


 

END MODULE m_eigen_diag
