!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_available_solvers
    IMPLICIT NONE
! list of constants indicating the different solvers
! the parameters are set to negative values to indicate that a particular solver is not compiled
! -solvers with numbers below 100 expect a non-distributed matrix
! -solvers with numbers 100-199 expect a distributed (scalapack-type) matrix
! -solvers with numbers higher than 200 can handle both
#ifdef CPP_ELPA
  INTEGER,PARAMETER:: diag_elpa=101
#else
  INTEGER,PARAMETER:: diag_elpa=-101
#endif
#ifdef CPP_ELEMENTAL
  INTEGER,PARAMETER:: diag_elemental=102
#else
  INTEGER,PARAMETER:: diag_elemental=-102
#endif
#ifdef CPP_SCALAPACK
  INTEGER,PARAMETER:: diag_scalapack=103
#else
  INTEGER,PARAMETER:: diag_scalapack=-103
#endif
#ifdef CPP_MAGMA
  INTEGER,PARAMETER:: diag_magma=6
#else
  INTEGER,PARAMETER:: diag_magma=-6
#endif
#ifdef CPP_CHASE
  INTEGER,PARAMETER:: diag_chase=207
#else
  INTEGER,PARAMETER:: diag_chase=-207
#endif
#ifdef CPP_CUSOLVER
  INTEGER,PARAMETER:: diag_cusolver=8
#else
  INTEGER,PARAMETER:: diag_cusolver=-8
#endif
  
  INTEGER,PARAMETER:: diag_lapack=4
#ifdef CPP_ELPA_ONENODE
  INTEGER,PARAMETER:: diag_elpa_1node=14
#else
  INTEGER,PARAMETER:: diag_elpa_1node=-14
#endif 
#ifdef CPP_SCALAPACK
  INTEGER,PARAMETER:: diag_debugout=201
#else
  INTEGER,PARAMETER:: diag_debugout=20
#endif  
  INTEGER,PARAMETER::diag_all_solver(9)=(/diag_elpa,diag_elemental,diag_scalapack,diag_magma,diag_chase,diag_cusolver,diag_lapack,diag_elpa_1node,diag_debugout/)
  
CONTAINS

  LOGICAL FUNCTION parallel_solver_available()
    parallel_solver_available=ANY(diag_all_solver>100)
  END FUNCTION parallel_solver_available

  FUNCTION select_solver(suggested_solver,parallel) RESULT(diag_solver)
    USE m_juDFT
    LOGICAL,INTENT(IN):: parallel
    INTEGER,INTENT(in):: suggested_solver
    INTEGER           :: diag_solver
    diag_solver=-99

    diag_solver=suggested_solver

    IF (suggested_solver==0) THEN
       !Determine the default solver
       IF (parallel) THEN
          diag_solver=MINVAL(diag_all_solver,mask=diag_all_solver>100)
       ELSE
          diag_solver=MINVAL(diag_all_solver,mask=diag_all_solver>0)
       ENDIF
       
       !check if a special solver was requested
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elpa")       diag_solver=diag_elpa
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elpa_1node") diag_solver=diag_elpa_1node
       IF (TRIM(juDFT_string_for_argument("-diag"))=="scalapack")  diag_solver=diag_scalapack
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elemental")  diag_solver=diag_elemental
       IF (TRIM(juDFT_string_for_argument("-diag"))=="lapack")     diag_solver=diag_lapack
       IF (TRIM(juDFT_string_for_argument("-diag"))=="magma")      diag_solver=diag_magma
       IF (TRIM(juDFT_string_for_argument("-diag"))=="chase")      diag_solver=diag_chase
       IF (TRIM(juDFT_string_for_argument("-diag"))=="cusolver")   diag_solver=diag_cusolver
       IF (TRIM(juDFT_string_for_argument("-diag"))=="debugout")   diag_solver=diag_debugout
       !Check if solver is possible
       IF (diag_solver<0)  CALL juDFT_error(&
            "You selected a solver for the eigenvalue problem that is not available",&
            hint="You most probably did not provide the appropriate libraries for compilation/linking")
       IF (diag_solver<100.AND.parallel) CALL judft_error(&
            "You selected an eigensolver that does not support distributed memory parallism",&
            hint="Try scalapack,elpa or another supported solver for parallel matrices")
       IF (diag_solver>100.AND.diag_solver<200.AND..NOT.parallel) CALL judft_error(&
            "You selected an eigensolver for matrices that are memory distributed",&
            hint="Try lapack, cusolver or another supported solver for non-distributed matrices")
    END IF

  END FUNCTION select_solver

END MODULE m_available_solvers
