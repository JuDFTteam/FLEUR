!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_matrix_pref
CONTAINS
   SUBROUTINE matrix_pref(fmpi, atoms, bmat, gvecPr, gvec, lapwPr, lapw, nvPr, nv, &
                          & iDtype, iDir, hmat_tmp, smat_tmp, hmat, smat, killcont)
      !> Decoratetes matrix elements of the form
      !! <\phi_{kG'q}|M|\phi_{kG}>
      !! with a prefactor i(G-G'-q).

      USE m_types

      IMPLICIT NONE

      TYPE(t_mpi),   INTENT(IN)    :: fmpi
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_lapw),  INTENT(IN)    :: lapwPr, lapw
      REAL,          INTENT(IN)    :: bmat(3, 3)
      INTEGER,       INTENT(IN)    :: gvecPr(:, :), gvec(:, :)
      INTEGER,       INTENT(IN)    :: nvPr, nv, iDtype, iDir, killcont(2)

      CLASS(t_mat),  INTENT(IN)    :: hmat_tmp, smat_tmp
      CLASS(t_mat),  INTENT(INOUT) :: hmat, smat

      INTEGER :: ikGPr, ikG, ikG0
      INTEGER :: iLo, l, imLo, ikLo, ikGLo, ikGLo0 ! LO stuff
      INTEGER :: iLoPr, lPr, imLoPr, ikLoPr, ikGLoPr ! LO stuff
      COMPLEX :: pref(3)

      !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
      !$OMP SHARED(fmpi, atoms, bmat, gvecPr, gvec, lapwPr, lapw, killcont) &
      !$OMP SHARED(nvPr, nv, iDir, hmat_tmp, smat_tmp, hmat, smat, iDtype) &
      !$OMP PRIVATE(ikGPr, ikG, ikG0, pref) &
      !$OMP PRIVATE(iLo, l, imLo, ikLo, ikGLo, ikGLo0) &
      !$OMP PRIVATE(iLoPr, lPr, imLoPr, ikLoPr, ikGLoPr)
      DO ikG = fmpi%n_rank + 1, nv + atoms%nlo(iDtype), fmpi%n_size
         ikG0 = (ikG-1) / fmpi%n_size + 1
         DO ikGPr = 1, nvPr + atoms%nlo(iDtype)
            IF (ikGPr<=nvPr.AND.ikG<=nv) THEN
               pref = gvec(:, ikG) + lapw%bkpt
               pref = pref - gvecPr(:, ikGPr) - lapwPr%bkpt
               pref = ImagUnit * MATMUL(pref, bmat)

               hmat%data_c(ikGPr, ikG0) = hmat%data_c(ikGPr, ikG0) &
                                      & + killcont(1) * pref(iDir) * hmat_tmp%data_c(ikGPr, ikG0)
               smat%data_c(ikGPr, ikG0) = smat%data_c(ikGPr, ikG0) &
                                      & + killcont(2) * pref(iDir) * smat_tmp%data_c(ikGPr, ikG0)
            ELSE IF (ikGPr<=nvPr.AND.ikG>nv) THEN
               iLo = ikG-nv
               l = atoms%llo(iLo, iDtype)
               DO imLo = 1, 2*l+1
                  ikLo = lapw%kvec(imLo,iLo,iDtype)
                  ikGLo = nv + lapw%index_lo(iLo,iDtype) + imLo
                  IF (MOD(ikGLo-1,fmpi%n_size) == fmpi%n_rank) THEN
                     ikGLo0 = (ikGLo-1)/fmpi%n_size+1

                     pref = gvec(:,ikLo) + lapw%bkpt
                     pref = pref - gvecPr(:, ikGPr) - lapwPr%bkpt
                     pref = ImagUnit * MATMUL(pref, bmat)

                     hmat%data_c(ikGPr, ikGLo0) = hmat%data_c(ikGPr, ikGLo0) &
                                            & + killcont(1) * pref(iDir) * hmat_tmp%data_c(ikGPr, ikGLo0)
                     smat%data_c(ikGPr, ikGLo0) = smat%data_c(ikGPr, ikGLo0) &
                                            & + killcont(2) * pref(iDir) * smat_tmp%data_c(ikGPr, ikGLo0)
                  END IF
               END DO
            ELSE IF (ikGPr>nvPr.AND.ikG<=nv) THEN
               iLoPr = ikGPr-nvPr
               lPr = atoms%llo(iLoPr, iDtype)
               DO imLoPr = 1, 2*lPr+1
                  ikLoPr = lapwPr%kvec(imLoPr,iLoPr,iDtype)
                  ikGLoPr = nvPr + lapwPr%index_lo(iLoPr,iDtype) + imLoPr

                  pref = gvec(:, ikG) + lapw%bkpt
                  pref = pref - gvec(:,ikLoPr) - lapwPr%bkpt
                  pref = ImagUnit * MATMUL(pref, bmat)

                  hmat%data_c(ikGLoPr, ikG0) = hmat%data_c(ikGLoPr, ikG0) &
                                           & + killcont(1) * pref(iDir) * hmat_tmp%data_c(ikGLoPr, ikG0)
                  smat%data_c(ikGLoPr, ikG0) = smat%data_c(ikGLoPr, ikG0) &
                                           & + killcont(2) * pref(iDir) * smat_tmp%data_c(ikGLoPr, ikG0)
               END DO
            ELSE
               iLoPr = ikGPr-nvPr
               lPr = atoms%llo(iLoPr, iDtype)
               iLo = ikG-nv
               l = atoms%llo(iLo, iDtype)
               DO imLo = 1, 2*l+1
                  ikLo = lapw%kvec(imLo,iLo,iDtype)
                  ikGLo = nv + lapw%index_lo(iLo,iDtype) + imLo
                  IF (MOD(ikGLo-1,fmpi%n_size) == fmpi%n_rank) THEN
                     ikGLo0 = (ikGLo-1)/fmpi%n_size+1

                     DO imLoPr = 1, 2*lPr+1
                        ikLoPr = lapwPr%kvec(imLoPr,iLoPr,iDtype)
                        ikGLoPr = nvPr + lapwPr%index_lo(iLoPr,iDtype) + imLoPr

                        pref = gvec(:,ikLo) + lapw%bkpt
                        pref = pref - gvec(:,ikLoPr) - lapwPr%bkpt
                        pref = ImagUnit * MATMUL(pref, bmat)

                        hmat%data_c(ikGLoPr, ikGLo0) = hmat%data_c(ikGLoPr, ikGLo0) &
                                                 & + killcont(1) * pref(iDir) * hmat_tmp%data_c(ikGLoPr, ikGLo0)
                        smat%data_c(ikGLoPr, ikGLo0) = smat%data_c(ikGLoPr, ikGLo0) &
                                                 & + killcont(2) * pref(iDir) * smat_tmp%data_c(ikGLoPr, ikGLo0)
                     END DO
                  END IF
               END DO
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE matrix_pref
END MODULE m_matrix_pref
