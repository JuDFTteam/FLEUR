!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_matrix_pref
CONTAINS
   SUBROUTINE matrix_pref(fmpi, bmat, gvecPr, gvec, kvecPr, kvec, nvPr, nv, &
                          & iDir, hmat_tmp, smat_tmp, hmat, smat)
      ! Calculates matrix elements of the form
      ! <\phi_{k'G'}|M|\phi_{kG}>
      ! for different use cases in the DFT/DFPT scf loop and operators M.
      !
      ! The spin loop and distinction is found on a higher level and translates
      ! into new switches.
      !
      ! DFT:
      ! M = \Theta_{IR} * (T + V) goes into hmat, M = \Theta_{IR} into smat.
      ! [vpw = \Theta_{IR} * V]
      ! [stars%ustep = \Theta_{IR}]
      ! [l_smat = F for offdiags]
      ! [iTkin = 0 for offdiags, 1 for l_useapw, 2 else]
      !
      ! DFPT:
      ! M = \Theta_{IR}^{(1)} * (T + V) + \Theta_{IR} * V^{(1)} goes into hmat,
      ! M = \Theta_{IR}^{(1)} into smat.
      ! [vpw = \Theta_{IR}^{(1)} * V + \Theta_{IR} * V^{(1)}]
      ! [stars%ustep = \Theta_{IR}^{(1)}]
      ! [l_smat = F for offdiags]
      ! [iTkin = 0 for offdiags, 1 else]

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
      REAL    :: pref(3)

      !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
      !$OMP SHARED(fmpi, bmat, gvecPr, gvec, kvecPr, kvec) &
      !$OMP SHARED(nvPr, nv, iDir, hmat_tmp, smat_tmp, hmat, smat) &
      !$OMP PRIVATE(ikGPr, ikG, ikG0, pref)
      DO ikG = fmpi%n_rank + 1, nv, fmpi%n_size
         ikG0 = (ikG-1) / fmpi%n_size + 1
         DO  ikGPr = 1, nvPr
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
