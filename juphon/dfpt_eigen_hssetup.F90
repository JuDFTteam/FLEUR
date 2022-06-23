!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_eigen_hssetup
CONTAINS
   SUBROUTINE dfpt_eigen_hssetup(isp, fmpi, fi, enpara, nococonv, starsq, &
                            ud, td, tdV1, vTot1, lapw, lapwq, iDir, iDtype, smat_final, hmat_final, nk, killcont)
      USE m_types
      USE m_types_mpimat
      USE m_types_gpumat
      USE m_dfpt_hs_int
      USE m_dfpt_hsmt
      USE m_dfpt_eigen_redist_matrix

      IMPLICIT NONE

      INTEGER,            INTENT(IN)     :: isp
      TYPE(t_mpi),        INTENT(IN)     :: fmpi
      type(t_fleurinput), INTENT(IN)     :: fi
      TYPE(t_stars),      INTENT(IN)     :: starsq
      TYPE(t_enpara),     INTENT(IN)     :: enpara
      TYPE(t_nococonv),   INTENT(IN)     :: nococonv
      TYPE(t_usdus),      INTENT(IN)     :: ud
      TYPE(t_tlmplm),     INTENT(IN)     :: td, tdV1
      TYPE(t_lapw),       INTENT(IN)     :: lapw, lapwq
      TYPE(t_potden),     INTENT(IN)     :: vTot1
      INTEGER,            INTENT(IN)     :: iDir, iDtype
      CLASS(t_mat), ALLOCATABLE, INTENT(INOUT)   :: smat_final, hmat_final
      INTEGER,      INTENT(IN)     :: nk, killcont(6)

      CLASS(t_mat), ALLOCATABLE :: smat(:, :), hmat(:, :)

      INTEGER :: i, j, nspins

      nspins = MERGE(2, 1, fi%noco%l_noco)
      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::smat(nspins, nspins), hmat(nspins, nspins))
      ELSE
         ALLOCATE (t_mpimat::smat(nspins, nspins), hmat(nspins, nspins))
      END IF

      DO i = 1, nspins
         DO j = 1, nspins
            CALL smat(i, j)%init(.FALSE., lapwq%nv(i) + fi%atoms%nlotot, lapw%nv(j) + fi%atoms%nlotot, fmpi%sub_comm, .false.)
            CALL hmat(i, j)%init(smat(i, j))
         END DO
      END DO

      CALL timestart("Interstitial part")
      CALL dfpt_hs_int(fi%noco, starsq, lapwq, lapw, fmpi, fi%cell%bbmat, isp, vTot1%pw_w, smat, hmat, killcont(1:3))
      CALL timestop("Interstitial part")

      CALL timestart("MT part")
      DO i = 1, nspins; DO j = 1, nspins
            !$acc enter data copyin(hmat(i,j),smat(i,j))
            !$acc enter data copyin(hmat(i,j)%data_r,smat(i,j)%data_r,hmat(i,j)%data_c,smat(i,j)%data_c)
      END DO; END DO
      CALL dfpt_hsmt(fi%atoms, fi%sym, enpara, isp, iDir, iDtype, fi%input, fmpi, fi%noco, nococonv, fi%cell, lapw, lapwq, ud, td, tdV1, hmat, smat, nk, killcont(4:6))
      DO i = 1, nspins; DO j = 1, nspins; if (hmat(1, 1)%l_real) THEN
            !$acc exit data copyout(hmat(i,j)%data_r,smat(i,j)%data_r) delete(hmat(i,j)%data_c,smat(i,j)%data_c)
            !$acc exist data delete(hmat(i,j),smat(i,j))
         ELSE
            !$acc exit data copyout(hmat(i,j)%data_c,smat(i,j)%data_c) delete(hmat(i,j)%data_r,smat(i,j)%data_r)
            !$acc exist data delete(hmat(i,j),smat(i,j))
         END IF; END DO; END DO
      CALL timestop("MT part")

      !Now copy the data into final matrix
      ! Collect the four fi%noco parts into a single matrix
      ! In collinear case only a copy is done
      ! In the parallel case also a redistribution happens
      ALLOCATE (smat_final, mold=smat(1, 1))
      ALLOCATE (hmat_final, mold=smat(1, 1))

      CALL timestart("Matrix redistribution")
      CALL dfpt_eigen_redist_matrix(fmpi, lapwq, lapw, fi%atoms, smat, smat_final)
      CALL dfpt_eigen_redist_matrix(fmpi, lapwq, lapw, fi%atoms, hmat, hmat_final, smat_final)
      CALL timestop("Matrix redistribution")

   END SUBROUTINE dfpt_eigen_hssetup
END MODULE m_dfpt_eigen_hssetup
