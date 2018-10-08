!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cusolver_diag
  USE m_types
  USE m_judft
#ifdef CPP_GPU  
  USE m_types_gpumat
#endif  
  IMPLICIT NONE
  PRIVATE
#ifdef CPP_GPU
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

#ifdef CPP_GPU
    INTEGER,PARAMETER:: max_sweeps=15
    REAL             :: tol=1E-7
    
    TYPE(t_gpumat)::smat_gpu,hmat_gpu
    
    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)
    SELECT TYPE(hmat)
    TYPE IS (t_mat)
       CALL hmat_gpu%init_from(hmat)
       CALL smat_gpu%init_from(smat)
       IF (hmat%l_real) THEN
          CALL cusolver_real(hmat_gpu%gpu_r,smat_gpu%gpu_r,smat%matsize1,ne,tol,max_sweeps,eig,z%data_r)
       ELSE
          CALL cusolver_complex(hmat_gpu%gpu_c,smat_gpu%gpu_c,smat%matsize1,ne,tol,max_sweeps,eig,z%data_c)
       END IF
    TYPE IS (t_gpumat)
       IF (hmat%l_real) THEN
          CALL cusolver_real(hmat%gpu_r,smat%gpu_r,smat%matsize1,ne,tol,max_sweeps,eig,z%data_r)
       ELSE
          CALL cusolver_complex(hmat%gpu_c,smat%gpu_c,smat%matsize1,ne,tol,max_sweeps,eig,z%data_c)
       END IF
       
    END SELECT
#endif
       
  END SUBROUTINE cusolver_diag

    
END MODULE m_cusolver_diag
