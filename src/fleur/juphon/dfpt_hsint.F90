!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_hs_int
CONTAINS
   ! Constructs the interstitial perturbed Hamiltonian and overlap matrix
   SUBROUTINE dfpt_hs_int(noco, juphon, starsq, lapwq, lapw, fmpi, bbmat, isp, vpw, hmat, smat, killcont)

      USE m_types
      USE m_hs_int_direct

      IMPLICIT NONE

      TYPE(t_noco),INTENT(IN)       :: noco
      TYPE(t_juphon),INTENT(IN)     :: juphon
      TYPE(t_stars),INTENT(IN)      :: starsq
      REAL, INTENT(IN)              :: bbmat(3, 3)
      TYPE(t_lapw),INTENT(IN)       :: lapwq, lapw
      TYPE(t_mpi),INTENT(IN)        :: fmpi
      INTEGER,INTENT(IN)            :: isp, killcont(3)
      COMPLEX,INTENT(IN)            :: vpw(:, :)
      CLASS(t_mat),INTENT(INOUT)    :: smat(:,:),hmat(:,:)

      INTEGER :: iSpinPr,iSpin, iMatPr, iMat, iTkin
      LOGICAL :: l_smat
      COMPLEX, ALLOCATABLE :: vpw_temp(:)

      IF (noco%l_noco.AND.isp.EQ.2) RETURN !was done already

      ALLOCATE(vpw_temp(SIZE(vpw, 1)))

      DO iSpinPr = MERGE(1, isp, noco%l_noco), MERGE(2, isp, noco%l_noco)
         ! co:
         ! iSpinPr = isp, SIZE(smat, 1) = 1 (?)
         ! noco:
         ! iSpinPr = 1...2, SIZE(smat, 1) = 2 (?)
         ! iispin = MIN(iSpinPr, SIZE(smat, 1))
         ! co:
         ! iispin = 1
         ! noco:
         ! iispin = 1...2
         ! --> alternative: iispin = MERGE(iSpinPr, 1, noco%l_noco) ?
         iMatPr = MERGE(iSpinPr, 1, noco%l_noco)
         DO iSpin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
            iMat = MERGE(iSpin, 1, noco%l_noco)
            iTkin = 0
            ! 1, 2, 3, 4 == 11, 22, 21, 12:
            IF ((iSpinPr.EQ.1).AND.(iSpin.EQ.1)) vpw_temp = vpw(:, 1) * killcont(1)
            IF ((iSpinPr.EQ.2).AND.(iSpin.EQ.2)) vpw_temp = vpw(:, 2) * killcont(1)
            IF ((iSpinPr.EQ.2).AND.(iSpin.EQ.1)) vpw_temp = vpw(:, 3) * killcont(1)
            IF ((iSpinPr.EQ.1).AND.(iSpin.EQ.2)) vpw_temp = vpw(:, 4) * killcont(1)

            l_smat = iSpinPr.EQ.iSpin
            IF (killcont(3)==0) l_smat = .FALSE.

            IF (iSpinPr.EQ.iSpin) iTkin = 2 * killcont(2)

            IF (.NOT.juphon%l_phonon) THEN
               l_smat = .FALSE.
               iTkin = 0
            END IF

            CALL hs_int_direct(fmpi, starsq, bbmat, lapwq%gvec(:, :, iSpinPr), lapw%gvec(:,:,iSpin), &
                             & lapwq%bkpt + lapwq%qphon, lapw%bkpt, lapwq%nv(iSpinPr), lapw%nv(iSpin), iTkin, 1, &
                             & l_smat, .TRUE., vpw_temp, hmat(iMatPr, iMat), smat(iMatPr, iMat))
         END DO
      END DO
   END SUBROUTINE dfpt_hs_int

   SUBROUTINE dfpt_dynmat_hs_int(noco, starsq, stars, lapwq, lapw, fmpi, bbmat, isp, theta1_pw0, theta1_pw, smat1, hmat1, smat1q, hmat1q, killcont)

      USE m_types
      USE m_hs_int_direct

      IMPLICIT NONE

      TYPE(t_noco),INTENT(IN)       :: noco
      TYPE(t_stars),INTENT(IN)      :: starsq, stars
      REAL, INTENT(IN)              :: bbmat(3, 3)
      TYPE(t_lapw),INTENT(IN)       :: lapwq, lapw
      TYPE(t_mpi),INTENT(IN)        :: fmpi
      INTEGER, INTENT(IN)           :: isp, killcont(4)
      COMPLEX, INTENT(IN)           :: theta1_pw0(:), theta1_pw(:)
      CLASS(t_mat),INTENT(INOUT)    :: smat1(:,:),hmat1(:,:),smat1q(:,:),hmat1q(:,:)!,smat2(:,:),hmat2(:,:)

      INTEGER :: iSpinPr,iSpin, iMatPr, iMat, iTkin
      LOGICAL :: l_smat
      COMPLEX, ALLOCATABLE :: vpw_temp(:), vpwq_temp(:)

      IF (noco%l_noco.AND.isp.EQ.2) RETURN !was done already

      ALLOCATE(vpw_temp(SIZE(stars%ustep, 1)))
      ALLOCATE(vpwq_temp(SIZE(starsq%ustep, 1)))

      DO iSpinPr = MERGE(1, isp, noco%l_noco), MERGE(2, isp, noco%l_noco)
         iMatPr = MERGE(iSpinPr, 1, noco%l_noco)
         DO iSpin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
            iMat = MERGE(iSpin, 1, noco%l_noco)
            iTkin = 0

            vpw_temp = CMPLX(0.0,0.0)
            vpwq_temp = CMPLX(0.0,0.0)

            l_smat = iSpinPr.EQ.iSpin
            IF (killcont(2)==0) l_smat = .FALSE.

            IF (iSpinPr.EQ.iSpin) iTkin = 2*killcont(1)

            CALL hs_int_direct(fmpi, stars, bbmat, lapw%gvec(:, :, iSpinPr), lapw%gvec(:,:,iSpin), &
                             & lapw%bkpt, lapw%bkpt, lapw%nv(iSpinPr), lapw%nv(iSpin), iTkin, 1, &
                             & l_smat, .TRUE., vpw_temp, hmat1(iMatPr, iMat), smat1(iMatPr, iMat), theta1_pw0)

            iTkin = 0
            l_smat = iSpinPr.EQ.iSpin
            IF (killcont(4)==0) l_smat = .FALSE.

            IF (iSpinPr.EQ.iSpin) iTkin = 2*killcont(3)
            CALL hs_int_direct(fmpi, starsq, bbmat, lapwq%gvec(:, :, iSpinPr), lapw%gvec(:,:,iSpin), &
                             & lapwq%bkpt + lapwq%qphon, lapw%bkpt, lapwq%nv(iSpinPr), lapw%nv(iSpin), iTkin, 1, &
                             & l_smat, .TRUE., vpwq_temp, hmat1q(iMatPr, iMat), smat1q(iMatPr, iMat), theta1_pw)
            ! Alternate form of the dynmat contribution:
            ! Calculate Theta2 as well
            !CALL hs_int_direct(fmpi, stars, bbmat, lapw%gvec(:, :, iSpinPr), lapw%gvec(:,:,iSpin), &
            !                 & lapw%bkpt, lapw%bkpt, lapw%nv(iSpinPr), lapw%nv(iSpin), iTkin, 1, &
            !                 & l_smat, .TRUE., vpw_temp, hmat2(iMatPr, iMat), smat2(iMatPr, iMat), theta1_pw0)
         END DO
      END DO
   END SUBROUTINE dfpt_dynmat_hs_int
END MODULE m_dfpt_hs_int
