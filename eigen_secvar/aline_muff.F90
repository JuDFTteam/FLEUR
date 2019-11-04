!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_alinemuff
  !************************************************************************
  !*                                                                      *
  !*     eigensystem-solver for moderatly-well converged potentials       *
  !*     a*z=e*b*z is transformed to h*z'=e*s*z' , whereby                *
  !*     h=C^T*a*C, s=C^T*b*C and z'=C^(-1)*z, when C is z of the last    *
  !*     iteration (lapw%nv*ne-array)                                          *
  !*     For ne<<lapw%nv the matrixsize is significantly reduced               *
  !*     aline uses ESSL-calls (use LAPACK's reduc3, tred3, bisect,       *
  !*     tinvit, trback and rebk3  if no ESSL available):                 *
  !*     SSPEV:  eigensystem-solver for symmetric, real packes h          *
  !*             here we have no s-matrix                                 *
  !*     For all eigenvalues are needed, SSPEV should perform better      *
  !*     then seclr4 (hope so)                                            *
  !*                                                     Gustav           *
  !*                                                                      *
  !************************************************************************
CONTAINS
  SUBROUTINE aline_muff(atoms,input,sym, cell, jsp,ne, usdus,td, bkpt,lapw, eig,z_r,z_c,realdata)

#include"cpp_double.h"

    USE m_hnonmuff
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_usdus),INTENT(IN)       :: usdus
    TYPE(t_lapw),INTENT(IN)        :: lapw
    TYPE(t_tlmplm),INTENT(IN)      :: td
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp,ne
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: bkpt(3)
    REAL,    INTENT (INOUT) :: eig(input%neig)

    REAL,    OPTIONAL,INTENT (INOUT) :: z_r(lapw%dim_nbasfcn(),ne)
    COMPLEX, OPTIONAL,INTENT (INOUT) :: z_c(lapw%dim_nbasfcn(),ne)
    LOGICAL,OPTIONAL,INTENT(IN):: realdata
    !     ..
    !     .. Local Scalars ..
    INTEGER i,info,j,ii
    !     ..
    !     .. Local Arrays ..
    REAL h(ne*(ne+1)/2),help(3*ne),z1(ne,ne)
    !     ..
    !     .. External Functions ..
    REAL CPP_BLAS_sdot
    EXTERNAL CPP_BLAS_sdot
    !     ..
    !     .. External Subroutines ..
    EXTERNAL CPP_LAPACK_ssygv
    LOGICAL l_real

    l_real=present(z_r)
    if (present(realdata)) l_real=realdata

    !     ..
    !---> initialize the hamiltonian and overlap matrix
       h = 0.0
       !---> add the diagonal (muffin-tin) terms
       DO i = 1,ne
          ii = (i-1)*i/2 + i
          h(ii) = eig(i)
       END DO

    !---> add the off-diagonal (non-muffin-tin) terms
    CALL h_nonmuff(atoms,input,sym, cell, jsp,ne, usdus,td, bkpt,lapw, h,l_real,z_r,z_c)

    !---> DIAGONALIZE THE HAMILTONIAN USING LIBRARY-ROUTINES
#ifdef CPP_ESSL
    !---> ESSL call, IBM AIX
    CALL CPP_LAPACK_sspev (21, h, eig,z1, ne,ne,help,3*ne)
#else
    !---> LAPACK call
    CALL CPP_LAPACK_sspev ('V','U',ne, h, eig,z1, ne,help, info)
    WRITE (6,FMT=8000) info
8000 FORMAT (' AFTER CPP_LAPACK_sspev: info=',i4)
#endif

    !---> store eigenvectors on array z
    DO i = 1,lapw%nv(jsp)
       if (l_real) THEN
          help(:ne)=z_r(i,:ne)
          DO j = 1,ne
             z_r(i,j) = CPP_BLAS_sdot(ne,help,1,z1(1,j),1)
          END DO
       else
          help(:ne)=z_c(i,:ne)
          DO j = 1,ne
             z_c(i,j) = CPP_BLAS_sdot(ne,help,1,z1(1,j),1)
          END DO
       endif
    END DO

  END SUBROUTINE aline_muff
END MODULE m_alinemuff
