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
! -solvers with numbers 100-199 can handle both distributed and non-distributed matrices
! -solvers with numbers higher than 200 expect a distributed (scalapack-type) matrix

!The default is choosen as: the minimum number >0 for non-distributed matrices
!                           the maximum number >100 for distributed matrices    
!1. solver
#ifdef CPP_CHASE
    INTEGER,PARAMETER:: diag_chase=207
#else
    INTEGER,PARAMETER:: diag_chase=-207
#endif
!2. solver
#ifdef CPP_ELEMENTAL
  INTEGER,PARAMETER:: diag_elemental=212
#else
  INTEGER,PARAMETER:: diag_elemental=-212
#endif
!3. solver
#ifdef CPP_SCALAPACK
  INTEGER,PARAMETER:: diag_scalapack=213
#else
  INTEGER,PARAMETER:: diag_scalapack=-213
#endif
!4. solver
#ifdef CPP_ELPA
  INTEGER,PARAMETER:: diag_elpa=214
#else
  INTEGER,PARAMETER:: diag_elpa=-214
#endif
!5.6. solver
#ifdef CPP_ELSI
INTEGER,PARAMETER:: diag_elsielpa=216
INTEGER,PARAMETER:: diag_elsichase=215
#else
  INTEGER,PARAMETER:: diag_elsielpa=-216
  INTEGER,PARAMETER:: diag_elsichase=-215
#endif
!7. solver
#ifdef CPP_MAGMA
  INTEGER,PARAMETER:: diag_magma=8
#else
  INTEGER,PARAMETER:: diag_magma=-8
#endif
!8. solver
  INTEGER,PARAMETER :: diag_stop=101 ! dummy solver that simply stops FLEUR
!9. solver
#ifdef CPP_CUSOLVER
INTEGER,PARAMETER:: diag_cusolver=7
#else
  INTEGER,PARAMETER:: diag_cusolver=-7
#endif
!10.11. solver
INTEGER,PARAMETER:: diag_lapack=10
INTEGER,PARAMETER:: diag_lapack_singlePrec=11
!12. solver  
#ifdef CPP_ELPA_ONENODE
  INTEGER,PARAMETER:: diag_elpa_1node=4
#else
  INTEGER,PARAMETER:: diag_elpa_1node=-4
#endif 
!13. solver
#ifdef CPP_SCALAPACK
  INTEGER,PARAMETER:: diag_debugout=120
#else
  INTEGER,PARAMETER:: diag_debugout=20
#endif  
!14.solver
INTEGER,PARAMETER:: diag_dummy=21


  INTEGER,PARAMETER::diag_all_solver(13)=(/diag_chase,diag_elemental,diag_scalapack,diag_elpa,diag_elsielpa,&
                     diag_elsichase,diag_magma,diag_stop,diag_cusolver,diag_lapack,diag_lapack_singlePrec,&
                     diag_elpa_1node,diag_debugout,diag_dummy/)
  
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
          diag_solver=MAXVAL(diag_all_solver,mask=diag_all_solver>100)
       ELSE
          diag_solver=MINVAL(diag_all_solver,mask=diag_all_solver<200.and.diag_all_solver>0)
       ENDIF
       
       !check if a special solver was requested
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elpa")       diag_solver=diag_elpa
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elpa_1node") diag_solver=diag_elpa_1node
       IF (TRIM(juDFT_string_for_argument("-diag"))=="scalapack")  diag_solver=diag_scalapack
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elemental")  diag_solver=diag_elemental
       IF (TRIM(juDFT_string_for_argument("-diag"))=="lapack")     diag_solver=diag_lapack
       IF (TRIM(juDFT_string_for_argument("-diag"))=="lapack_singlePrec")     diag_solver=diag_lapack_singlePrec
       IF (TRIM(juDFT_string_for_argument("-diag"))=="magma")      diag_solver=diag_magma
       IF (TRIM(juDFT_string_for_argument("-diag"))=="chase")      diag_solver=diag_chase
       IF (TRIM(juDFT_string_for_argument("-diag"))=="cusolver")   diag_solver=diag_cusolver
       IF (TRIM(juDFT_string_for_argument("-diag"))=="debugout")   diag_solver=diag_debugout
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elsielpa")   diag_solver=diag_elsielpa
       IF (TRIM(juDFT_string_for_argument("-diag"))=="elsichase")  diag_solver=diag_elsichase
       IF (TRIM(juDFT_string_for_argument("-diag"))=="dummy")      diag_solver=diag_dummy
       IF (TRIM(juDFT_string_for_argument("-diag"))=="stop")       diag_solver=diag_stop
       
       !Check if solver is possible
       IF (diag_solver<0)  CALL juDFT_error(&
            "You selected a solver for the eigenvalue problem that is not available",&
            hint="You most probably did not provide the appropriate libraries for compilation/linking")
       IF (diag_solver<100.AND.parallel) CALL judft_error(&
            "You selected an eigensolver that does not support distributed memory parallism",&
            hint="Try scalapack,elpa or another supported solver for parallel matrices")
       IF (diag_solver>200.AND..NOT.parallel) CALL judft_error(&
            "You selected an eigensolver for matrices that are memory distributed",&
            hint="Try lapack, cusolver or another supported solver for non-distributed matrices")
    END IF

  END FUNCTION select_solver

END MODULE m_available_solvers
