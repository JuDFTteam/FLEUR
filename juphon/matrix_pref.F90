!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_matrix_pref
CONTAINS
   SUBROUTINE matrix_pref(fmpi, bmat, gvecPr, gvec, kvecPr, kvec, nvPr, nv, &
                          & iDir, hmat_tmp, smat_tmp, hmat, smat)
      !> Decoratetes matrix elements of the form
      !! <\phi_{kG'q}|M|\phi_{kG}>
      !! with a prefactor i(G-G'-q).

      USE m_types

      IMPLICIT NONE

      TYPE(t_mpi),   INTENT(IN)    :: fmpi
      REAL,          INTENT(IN)    :: bmat(3, 3)
      INTEGER,       INTENT(IN)    :: gvecPr(:, :), gvec(:, :)
      REAL,          INTENT(IN)    :: kvecPr(3), kvec(3)
      INTEGER,       INTENT(IN)    :: nvPr, nv, iDir

      CLASS(t_mat),  INTENT(IN)    :: hmat_tmp, smat_tmp
      CLASS(t_mat),  INTENT(INOUT) :: hmat, smat

      INTEGER :: ikGPr, ikG, ikG0
      COMPLEX :: pref(3)

      !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
      !$OMP SHARED(fmpi, bmat, gvecPr, gvec, kvecPr, kvec) &
      !$OMP SHARED(nvPr, nv, iDir, hmat_tmp, smat_tmp, hmat, smat) &
      !$OMP PRIVATE(ikGPr, ikG, ikG0, pref)
      DO ikG = fmpi%n_rank + 1, nv, fmpi%n_size
         ikG0 = (ikG-1) / fmpi%n_size + 1
         DO ikGPr = 1, nvPr
            pref = gvec(:, ikG) + kvec
            pref = pref - gvecPr(:, ikGPr) - kvecPr
            pref = ImagUnit * MATMUL(pref, bmat)

            hmat%data_c(ikGPr, ikG0) = hmat%data_c(ikGPr, ikG0) &
                                   & + pref(iDir) * hmat_tmp%data_c(ikGPr, ikG0)
            smat%data_c(ikGPr, ikG0) = smat%data_c(ikGPr, ikG0) &
                                   & + pref(iDir) * smat_tmp%data_c(ikGPr, ikG0)
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE matrix_pref
END MODULE m_matrix_pref
