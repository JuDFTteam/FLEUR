!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cusolver_diag
  USE m_types_mat
  USE m_types_mpimat
  USE m_judft
#ifdef CPP_GPU  
  USE m_types_gpumat
#endif  
  IMPLICIT NONE
  PRIVATE
#ifdef CPP_CUSOLVER
  INTERFACE
     SUBROUTINE cusolver_real(H,S,n,ne,tol,max_sweeps,eig,z) BIND(C,name="cusolver_real") 
      USE iso_c_binding
      IMPLICIT NONE
      REAL(c_double)         :: H(*),S(*)
      INTEGER(c_int),VALUE   :: n,ne,max_sweeps
      REAL(c_double),VALUE   :: tol
      REAL(c_double)         :: eig(*),z(*)
    END SUBROUTINE cusolver_real
 END INTERFACE
 INTERFACE
    SUBROUTINE cusolver_complex(H,S,n,ne,tol,max_sweeps,eig,z) BIND(C,name="cusolver_real") 
      USE iso_c_binding
      IMPLICIT NONE
      COMPLEX(c_double)      :: H(*),S(*)
      INTEGER(c_int),VALUE   :: n,ne,max_sweeps
      REAL(c_double),VALUE   :: tol
      REAL(c_double)         :: eig(*)
      COMPLEX(c_double)      :: z(*)
    END SUBROUTINE cusolver_complex
 END INTERFACE
#endif
 PUBLIC cusolver_diag

CONTAINS
  SUBROUTINE cusolver_diag(hmat,smat,ne,eig,zmat)
    !Simple driver to solve Generalized Eigenvalue Problem using CuSolverDN
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT) :: hmat,smat
    INTEGER,INTENT(INOUT)      :: ne
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)    :: zmat
    REAL,INTENT(OUT)           :: eig(:)

#ifdef CPP_CUSOLVER
    INTEGER,PARAMETER:: max_sweeps=15
    REAL             :: tol=1E-7
    
    
    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)
    IF (hmat%l_real) THEN
       CALL cusolver_real(hmat%data_r,smat%data_r,smat%matsize1,ne,tol,max_sweeps,eig,zmat%data_r)
    ELSE
       CALL cusolver_complex(hmat%data_c,smat%data_c,smat%matsize1,ne,tol,max_sweeps,eig,zmat%data_c)
    END IF
#endif
       
  END SUBROUTINE cusolver_diag

    
END MODULE m_cusolver_diag
