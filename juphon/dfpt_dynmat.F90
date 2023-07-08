!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_dynmat
   USE m_types
   USE m_constants

IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt_dynmat_row(fi, stars, starsq, sphhar, xcpot, nococonv, hybdat, fmpi, qpts, iQ, iDtype_row, iDir_row, &
                              eig_id, dfpt_eig_id, dfpt_eig_id2, enpara, mpdata, results, results1, l_real, &
                              rho, vTot, grRho3, grVext3, grVC3, grVtot3, denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, dyn_row, &
                              l_dynMat0, l_dynMatq, E2ndOrdII, q_eig_id)
      USE m_step_function
      USE m_convol
      USE m_dfpt_vgen
      USE m_npy

      TYPE(t_fleurinput), INTENT(IN)    :: fi
      TYPE(t_stars),      INTENT(IN)    :: stars, starsq
      TYPE(t_sphhar),     INTENT(IN)    :: sphhar
      CLASS(t_xcpot),     INTENT(IN)    :: xcpot
      TYPE(t_nococonv),   INTENT(IN)    :: nococonv
      TYPE(t_hybdat),     INTENT(INOUT) :: hybdat
      TYPE(t_mpi),        INTENT(IN)    :: fmpi
      TYPE(t_kpts),       INTENT(IN)    :: qpts

      TYPE(t_potden), INTENT(INOUT) :: rho, vTot
      TYPE(t_potden), INTENT(INOUT) :: grRho3(3), grVext3(3), grVC3(3), grVtot3(3)
      TYPE(t_potden), INTENT(INOUT) :: denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im

      TYPE(t_enpara),   INTENT(INOUT) :: enpara
      TYPE(t_mpdata),   INTENT(INOUT) :: mpdata
      TYPE(t_results),  INTENT(INOUT) :: results, results1

      LOGICAL, INTENT(IN) :: l_real, l_dynMat0, l_dynMatq

      INTEGER, INTENT(IN) :: iQ, iDtype_row, iDir_row, eig_id, dfpt_eig_id, dfpt_eig_id2

      COMPLEX, INTENT(INOUT) :: dyn_row(:)

      COMPLEX, INTENT(IN) :: E2ndOrdII(:,:)
     
      INTEGER, OPTIONAL, INTENT(IN) :: q_eig_id

      TYPE(t_fftgrid) :: fftgrid_dummy
      TYPE(t_potden)  :: rho_dummy, rho1_dummy, vExt1, vExt1Im
      TYPE(t_hub1data) :: hub1data

      INTEGER :: col_index, row_index, iDtype_col, iDir_col, iType, iDir, iSpin, ierr
      COMPLEX :: tempval
      LOGICAL :: bare_mode

      REAL :: qvec(3)

      COMPLEX, ALLOCATABLE :: dyn_row_HF(:), dyn_row_eigen(:), dyn_row_int(:)
      COMPLEX, ALLOCATABLE :: theta1full(:, :, :), theta1full0(:, :, :)!, theta2(:, :, :)
      COMPLEX, ALLOCATABLE :: theta1_pw(:, :, :), theta1_pw0(:, :, :)!,theta2_pw(:, :, :)
      COMPLEX, ALLOCATABLE :: pww(:), pwwq(:), pww2(:), pwwq2(:)
      COMPLEX, ALLOCATABLE :: rho_pw(:), denIn1_pw(:)
      REAL,    ALLOCATABLE :: rho_mt(:,:,:), grRho_mt(:,:,:), denIn1_mt(:,:,:), denIn1_mt_Im(:,:,:)

      type(t_fft) :: fft

      bare_mode = .FALSE.

      ALLOCATE(dyn_row_HF(SIZE(dyn_row)), dyn_row_eigen(SIZE(dyn_row)), dyn_row_int(SIZE(dyn_row)))
      ALLOCATE(theta1full(0:27*starsq%mx1*starsq%mx2*starsq%mx3-1,fi%atoms%ntype,3))
      ALLOCATE(theta1full0(0:27*stars%mx1*stars%mx2*stars%mx3-1,fi%atoms%ntype,3))
      ALLOCATE(theta1_pw(starsq%ng3,fi%atoms%ntype,3))
      ALLOCATE(theta1_pw0(stars%ng3,fi%atoms%ntype,3))
      ALLOCATE(pww(stars%ng3),pwwq(starsq%ng3))
      ALLOCATE(pww2(stars%ng3),pwwq2(starsq%ng3))

      ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
      ALLOCATE(denIn1_pw(starsq%ng3),rho_pw(stars%ng3))
      ALLOCATE(rho_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),grRho_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))

      CALL rho_dummy%copyPotDen(rho)
      CALL rho1_dummy%copyPotDen(denIn1)
      CALL rho_dummy%resetPotDen()
      CALL rho1_dummy%resetPotDen()

      qvec = qpts%bk(:, iQ)

      theta1full  = CMPLX(0.0,0.0)
      theta1full0 = CMPLX(0.0,0.0)

      CALL stepf_analytical(fi%sym, starsq, fi%atoms, fi%input, fi%cell, fmpi, fftgrid_dummy, qvec, iDtype_row, iDir_row, 1, theta1full)
      CALL stepf_analytical(fi%sym, stars, fi%atoms, fi%input, fi%cell, fmpi, fftgrid_dummy, [0.0,0.0,0.0], iDtype_row, iDir_row, 1, theta1full0)
      !CALL stepf_analytical(fi%sym, stars, fi%atoms, fi%input, fi%cell, fmpi, fftgrid_dummy, [0.0,0.0,0.0], iDtype_row, iDir_row, 2, theta2)

      DO iType = 1, fi%atoms%ntype
         DO iDir = 1, 3
            fftgrid_dummy%grid = theta1full(0:, iType, iDir)
            CALL fftgrid_dummy%takeFieldFromGrid(starsq, theta1_pw(:, iType, iDir))
            theta1_pw(:, iType, iDir) = theta1_pw(:, iType, iDir) * 3 * starsq%mx1 * 3 * starsq%mx2 * 3 * starsq%mx3
            CALL fftgrid_dummy%perform_fft(forward=.false.)
            theta1full(0:, iType, iDir) = fftgrid_dummy%grid

            fftgrid_dummy%grid = theta1full0(0:, iType, iDir)
            CALL fftgrid_dummy%takeFieldFromGrid(stars, theta1_pw0(:, iType, iDir))
            theta1_pw0(:, iType, iDir) = theta1_pw0(:, iType, iDir) * 3 * stars%mx1 * 3 * stars%mx2 * 3 * stars%mx3
            CALL fftgrid_dummy%perform_fft(forward=.false.)
            theta1full0(0:, iType, iDir) = fftgrid_dummy%grid
         END DO
      END DO

      row_index = 3 * (iDtype_row - 1) + iDir_row

      dyn_row       = CMPLX(0.0,0.0)
      dyn_row_HF    = CMPLX(0.0,0.0)
      dyn_row_eigen = CMPLX(0.0,0.0)
      dyn_row_int   = CMPLX(0.0,0.0)

      denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
      denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
      ! Get "full" denIn1:
      denIn1_mt(:,0:,iDtype_row) = denIn1_mt(:,0:,iDtype_row) - (grRho3(iDir_row)%mt(:,0:,iDtype_row,1)+grRho3(iDir_row)%mt(:,0:,iDtype_row,fi%input%jspins))/(3.0-fi%input%jspins)
      denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)

      DO iDtype_col = 1, fi%atoms%ntype
         DO iDir_col = 1, 3
            IF (fmpi%irank==0) write(9989,*) "------------------"
            IF (fmpi%irank==0) write(9989,*) "Atom:", iDtype_col, "Direction", iDir_col
            tempval = CMPLX(0.0,0.0)
            col_index = 3 * (iDtype_col - 1) + iDir_col

            ! First calculate the HF contributions.
            ! \rho(1)V_{ext}(1) integral over whole unit cell

            ! Get V_{ext}(1) for \alpha, i with gradient cancellation
            CALL vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
            CALL vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
            CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_dummy,vTot,&
                           starsq,rho1_dummy,vExt1,.FALSE.,vExt1Im,rho1_dummy,iDtype_col,iDir_col,[0,0])

            ! IR integral:
            pwwq2 = CMPLX(0.0,0.0)
            CALL dfpt_convol_big(1, starsq, stars, vExt1%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pwwq2)
            CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval)
            dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
            IF (fmpi%irank==0) write(9989,*) "IR rho1 V1ext new             ", tempval
            tempval = CMPLX(0.0,0.0)

            ! MT integral:
            ! If we use gradient cancellation, remove it from rho1
            IF (.NOT.bare_mode) denIn1_mt(:,0:,iDtype_row) = &
                                denIn1_mt(:,0:,iDtype_row) + &
                                (grRho3(iDir_row)%mt(:,0:,iDtype_row,1)+grRho3(iDir_row)%mt(:,0:,iDtype_row,fi%input%jspins))/(3.0-fi%input%jspins)
            DO iType = 1, fi%atoms%ntype
               IF (fmpi%irank==0) write(9989,*) "Loop atom:", iType
               CALL dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval)
               dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
               IF (fmpi%irank==0) write(9989,*) "    MT rho1 V1ext                 ", tempval
               tempval = CMPLX(0.0,0.0)
            END DO
            IF (fmpi%irank==0) write(9989,*) "End atom loop"
            IF (.NOT.bare_mode) denIn1_mt(:,0:,iDtype_row) = &
                                denIn1_mt(:,0:,iDtype_row) - &
                                (grRho3(iDir_row)%mt(:,0:,iDtype_row,1)+grRho3(iDir_row)%mt(:,0:,iDtype_row,fi%input%jspins))/(3.0-fi%input%jspins)

            ! Various V_ext integrals:
            ! IR:
            rho_pw = (rho%pw(:,1)+rho%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
            pwwq2 = CMPLX(0.0,0.0)
            CALL dfpt_convol_big(2, stars, starsq, rho_pw, theta1full(0:, iDtype_row, iDir_row), pwwq2)
            CALL dfpt_int_pw(starsq, fi%cell, pwwq2, vExt1%pw(:,1), tempval)
            dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
            IF (fmpi%irank==0) write(9989,*) "IR Theta1 rho V1ext new       ", tempval
            tempval = CMPLX(0.0,0.0)

            ! MT:
            ! TODO: Apply baremode correction to MT correction and corresponding 2nd integral. They
            !       are the last two *big* terms. 
            IF (.NOT.bare_mode) denIn1_mt(:,0:,iDtype_row) = &
                                denIn1_mt(:,0:,iDtype_row) + &
                                (grRho3(iDir_row)%mt(:,0:,iDtype_row,1)+grRho3(iDir_row)%mt(:,0:,iDtype_row,fi%input%jspins))/(3.0-fi%input%jspins)
            grRho_mt = grVC3(iDir_col)%mt(:,0:,:,1)
            CALL dfpt_int_mt(fi%atoms, sphhar, fi%sym, iDtype_col, denIn1_mt, denIn1_mt_Im, grRho_mt, 0*grRho_mt, tempval)
            dyn_row_int(col_index) = dyn_row_int(col_index) + tempval
            IF (fmpi%irank==0) write(9989,*) "MT rho1 grVC                  ", tempval
            tempval = CMPLX(0.0,0.0)
            IF (.NOT.bare_mode) denIn1_mt(:,0:,iDtype_row) = &
                                denIn1_mt(:,0:,iDtype_row) - &
                                (grRho3(iDir_row)%mt(:,0:,iDtype_row,1)+grRho3(iDir_row)%mt(:,0:,iDtype_row,fi%input%jspins))/(3.0-fi%input%jspins)

            IF (.NOT.bare_mode) vC1%mt(:,0:,iDtype_row,1) = &
                                vC1%mt(:,0:,iDtype_row,1) + &
                                grVC3(iDir_row)%mt(:,0:,iDtype_row,1)
            grRho_mt = -(grRho3(iDir_col)%mt(:,0:,:,1)+grRho3(iDir_col)%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            CALL dfpt_int_mt(fi%atoms, sphhar, fi%sym, iDtype_col, vC1%mt(:,0:,:,1), vC1Im%mt(:,0:,:,1), grRho_mt, 0*grRho_mt, tempval)
            dyn_row_int(col_index) = dyn_row_int(col_index) + tempval
            IF (fmpi%irank==0) write(9989,*) "MT grRho V1C                  ", tempval
            tempval = CMPLX(0.0,0.0)
            IF (.NOT.bare_mode) vC1%mt(:,0:,iDtype_row,1) = &
                                vC1%mt(:,0:,iDtype_row,1) - &
                                grVC3(iDir_row)%mt(:,0:,iDtype_row,1)

            IF (.NOT.bare_mode) THEN
               IF (.NOT.bare_mode) vExt1%mt(:,0:,iDtype_col,:) = vExt1%mt(:,0:,iDtype_col,:) + grVext3(iDir_col)%mt(:,0:,iDtype_col,:)
               grRho_mt = -(grRho3(iDir_row)%mt(:,0:,:,1)+grRho3(iDir_row)%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
               CALL dfpt_int_mt(fi%atoms, sphhar, fi%sym, iDtype_row, grRho_mt, 0*grRho_mt, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval)
               dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
               IF (fmpi%irank==0) write(9989,*) "MT correction grRho V1ext     ", tempval
               tempval = CMPLX(0.0,0.0)
               IF (.NOT.bare_mode) vExt1%mt(:,0:,iDtype_col,:) = vExt1%mt(:,0:,iDtype_col,:) - grVext3(iDir_col)%mt(:,0:,iDtype_col,:)
            END IF

            ! SF:
            rho_mt = (rho%mt(:,0:,:,1)+rho%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            CALL dfpt_int_mt_sf(fi%atoms, sphhar, fi%sym, iDir_row, iDtype_row, rho_mt, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval)
            dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
            IF (fmpi%irank==0) write(9989,*) "SF rho Vext1                  ", tempval
            tempval = CMPLX(0.0,0.0)

            IF (.NOT.bare_mode) vC1%mt(:,0:,iDtype_row,1) = &
                                vC1%mt(:,0:,iDtype_row,1) + &
                                grVC3(iDir_row)%mt(:,0:,iDtype_row,1)
            CALL dfpt_int_mt_sf(fi%atoms, sphhar, fi%sym, iDir_col, iDtype_col, rho_mt, vC1%mt(:,0:,:,1), -vC1Im%mt(:,0:,:,1), tempval)
            dyn_row_int(col_index) = dyn_row_int(col_index) + tempval
            IF (fmpi%irank==0) write(9989,*) "SF rho VC1                    ", tempval
            tempval = CMPLX(0.0,0.0)
            IF (.NOT.bare_mode) vC1%mt(:,0:,iDtype_row,1) = &
                                vC1%mt(:,0:,iDtype_row,1) - &
                                grVC3(iDir_row)%mt(:,0:,iDtype_row,1)

            ! Miscellaneous integrals:
            rho_pw = (rho%pw(:,1)+rho%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
            pwwq2 = CMPLX(0.0,0.0)
            CALL dfpt_convol_big(2, stars, starsq, rho_pw, theta1full(0:, iDtype_col, iDir_col), pwwq2)
            CALL dfpt_int_pw(starsq, fi%cell, vC1%pw(:,1), pwwq2, tempval)
            dyn_row_int(col_index) = dyn_row_int(col_index) + tempval
            IF (fmpi%irank==0) write(9989,*) "IR V1C rho Theta1 new         ", tempval
            tempval = CMPLX(0.0,0.0)

            DO iSpin = 1, fi%input%jspins
               IF (fmpi%irank==0) write(9989,*) "Loop spin:", iSpin
               ! TODO: Ensure, that vTot/denIn1 is diagonal here, not 2x2.
               pwwq2 = CMPLX(0.0,0.0)
               CALL dfpt_convol_big(2, stars, starsq, vTot%pw(:, iSpin), theta1full(0:, iDtype_col, iDir_col), pwwq2)
               CALL dfpt_int_pw(starsq, fi%cell, denIn1%pw(:,iSpin), pwwq2, tempval)
               dyn_row_int(col_index) = dyn_row_int(col_index) + tempval
               IF (fmpi%irank==0) write(9989,*) "    IR rho1 vTot Theta1 new       ", tempval
               tempval = CMPLX(0.0,0.0)
            END DO
            IF (fmpi%irank==0) write(9989,*) "End spin loop"

            IF (iDtype_row==iDtype_col) THEN
               CALL vExt1%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
               CALL vExt1Im%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
               ! Get V_{ext}(1) for \alpha, i, q=0 with gradient cancellation
               CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                              fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_dummy,vTot,&
                              stars,rho_dummy,vExt1,.FALSE.,vExt1Im,rho_dummy,iDtype_col,iDir_col,[0,0])

               ! Integrals:
               rho_pw = (grRho3(iDir_row)%pw(:,1)+grRho3(iDir_row)%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
               pww2 = CMPLX(0.0,0.0)
               CALL dfpt_convol_big(1, stars, stars, vExt1%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pww2)
               CALL dfpt_int_pw(stars, fi%cell, rho_pw, pww2, tempval)
               dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
               IF (fmpi%irank==0) write(9989,*) "IR grRho V1ext0 new           ", tempval
               tempval = CMPLX(0.0,0.0)
               rho_mt = (rho%mt(:,0:,:,1)+rho%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
               rho_pw = -(rho%pw(:,1)+rho%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
               DO iType = 1, fi%atoms%ntype
                  IF (fmpi%irank==0) write(9989,*) "Loop atom:", iType
                  pww2 = CMPLX(0.0,0.0)
                  CALL dfpt_convol_big(1, stars, stars, rho_pw, theta1full0(0:, iType, iDir_row), pww2)
                  CALL dfpt_int_pw(stars, fi%cell, pww2, vExt1%pw(:,1), tempval)
                  dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
                  IF (fmpi%irank==0) write(9989,*) "   IR Theta1 rho V1ext0 new      ", tempval
                  tempval = CMPLX(0.0,0.0)

                  IF (.NOT.bare_mode) vExt1%mt(:,0:,iDtype_col,:) = vExt1%mt(:,0:,iDtype_col,:) + grVext3(iDir_col)%mt(:,0:,iDtype_col,:)
                  grRho_mt = (grRho3(iDir_row)%mt(:,0:,:,1)+grRho3(iDir_row)%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
                  CALL dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, grRho_mt, 0*grRho_mt, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval)
                  dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
                  IF (fmpi%irank==0) write(9989,*) "    MT grRho V1ext0               ", tempval
                  tempval = CMPLX(0.0,0.0)
                  IF (.NOT.bare_mode) vExt1%mt(:,0:,iDtype_col,:) = vExt1%mt(:,0:,iDtype_col,:) - grVext3(iDir_col)%mt(:,0:,iDtype_col,:)

                  CALL dfpt_int_mt_sf(fi%atoms, sphhar, fi%sym, iDir_row, iType, -rho_mt, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval)
                  dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
                  IF (fmpi%irank==0) write(9989,*) "    SF rho V1ext0                 ", tempval
                  tempval = CMPLX(0.0,0.0)
               END DO
               IF (fmpi%irank==0) write(9989,*) "End atom loop"

               rho_pw = (rho%pw(:,1)+rho%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
               pww2 = CMPLX(0.0,0.0)
               CALL dfpt_convol_big(1, stars, stars, rho_pw, theta1full0(0:, iDtype_col, iDir_col), pww2)
               CALL dfpt_int_pw(stars, fi%cell, grVC3(iDir_row)%pw(:,1), pww2, tempval)
               dyn_row_int(col_index) = dyn_row_int(col_index) + tempval
               IF (fmpi%irank==0) write(9989,*) "IR grVC rho Theta1 new        ", tempval
               tempval = CMPLX(0.0,0.0)

               DO iSpin = 1, fi%input%jspins
                  IF (fmpi%irank==0) write(9989,*) "Loop spin:", iSpin
                  ! TODO: Ensure, that vTot/gradrho is diagonal here, not 2x2 [NOCO].
                  pww2 = CMPLX(0.0,0.0)
                  CALL dfpt_convol_big(1, stars, stars, vTot%pw(:,iSpin), theta1full0(0:,iDtype_col,iDir_col), pww2)
                  CALL dfpt_int_pw(stars, fi%cell, grRho3(iDir_row)%pw(:,iSpin), pww2, tempval)
                  dyn_row_int(col_index) = dyn_row_int(col_index) + tempval
                  IF (fmpi%irank==0) write(9989,*) "    IR grRho vTot Theta1 new      ", tempval
                  tempval = CMPLX(0.0,0.0)
               END DO
               IF (fmpi%irank==0) write(9989,*) "End spin loop"
            END IF

            IF (fmpi%irank==0) write(9990,*) qvec, iDtype_row, iDir_row, iDtype_col, iDir_col
            IF (fmpi%irank==0) write(9990,*) "HF   :", conjg(dyn_row_HF(col_index))
            IF (fmpi%irank==0) write(9990,*) "int  :", conjg(dyn_row_int(col_index))

            ! Calculate the contributions to the dynamical matrix that stem
            ! from terms related to occupation numbers and the eigenenergies.
            IF (.NOT.PRESENT(q_eig_id)) THEN
               CALL dfpt_dynmat_eigen(fi, results, results1, xcpot, fmpi, mpdata, hybdat, enpara, nococonv, &
                                      stars, starsq, sphhar, rho, hub1data, vTot, vTot, vTot1, vTot1Im, &
                                      eig_id, dfpt_eig_id, dfpt_eig_id2, iDir_col, iDtype_col, iDir_row, iDtype_row, &
                                      theta1_pw0(:,iDtype_col,iDir_col), theta1_pw(:,iDtype_col,iDir_col), &
                                      qvec, l_real, dyn_row_eigen(col_index),[1,1,1,1,1,1,1,1,1,1,1])
            ELSE
               CALL dfpt_dynmat_eigen(fi, results, results1, xcpot, fmpi, mpdata, hybdat, enpara, nococonv, &
                                   stars, starsq, sphhar, rho, hub1data, vTot, vTot, vTot1, vTot1Im, &
                                   eig_id, dfpt_eig_id, dfpt_eig_id2, iDir_col, iDtype_col, iDir_row, iDtype_row, &
                                   theta1_pw0(:,iDtype_col,iDir_col), theta1_pw(:,iDtype_col,iDir_col), &
                                   qvec, l_real, dyn_row_eigen(col_index),[1,1,1,1,1,1,1,1,1,1,1],q_eig_id) 
            END IF
            IF (fmpi%irank==0) write(9990,*) "eigen:", dyn_row_eigen(col_index)
            IF (fmpi%irank==0) write(9989,*) "Eii2:", conjg(E2ndOrdII(3*(iDtype_row-1)+iDir_row,3*(iDtype_col-1)+iDir_col))
            IF (fmpi%irank==0) write(9990,*) "Eii2:", conjg(E2ndOrdII(3*(iDtype_row-1)+iDir_row,3*(iDtype_col-1)+iDir_col))
         END DO
      END DO

      dyn_row = conjg(dyn_row_HF) + conjg(dyn_row_int) + dyn_row_eigen

   END SUBROUTINE dfpt_dynmat_row

   SUBROUTINE dfpt_int_pw(stars, cell, pw_conj, pw_pure, pw_int)
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_cell),  INTENT(IN) :: cell

      COMPLEX, INTENT(IN)  :: pw_conj(:), pw_pure(:)

      COMPLEX, INTENT(INOUT) :: pw_int

      pw_int = pw_int + cell%omtil * DOT_PRODUCT(pw_conj(:stars%ng3),pw_pure(:stars%ng3))
   END SUBROUTINE dfpt_int_pw

   SUBROUTINE dfpt_int_mt(atoms, sphhar, sym, nat, mt_conj, mt_conj_im, mt_pure, mt_pure_im, mt_int)
      USE m_intgr, ONLY: intgr3, intgr3LinIntp

      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_sphhar), INTENT(IN) :: sphhar
      TYPE(t_sym),    INTENT(IN) :: sym

      INTEGER, INTENT(IN)  :: nat

      REAL, INTENT(IN)  :: mt_conj(:,0:,:), mt_conj_im(:,0:,:), mt_pure(:,0:,:), mt_pure_im(:,0:,:)

      COMPLEX, INTENT(INOUT) :: mt_int

      REAL    :: dpdot_re, dpdot_im
      COMPLEX :: tmt
      INTEGER :: j, lh

      REAL :: dpj_re(atoms%jmtd), dpj_im(atoms%jmtd)

      tmt = CMPLX(0.0,0.0)
      DO lh = 0, sphhar%nlh(sym%ntypsy(nat))
         DO j = 1, atoms%jri(nat)
            dpj_re(j) = mt_conj(j,lh,nat)*mt_pure(j,lh,nat)+mt_conj_im(j,lh,nat)*mt_pure_im(j,lh,nat)
            dpj_im(j) = mt_conj(j,lh,nat)*mt_pure_im(j,lh,nat)-mt_conj_im(j,lh,nat)*mt_pure(j,lh,nat)
         END DO
         CALL intgr3(dpj_re,atoms%rmsh(1,nat),atoms%dx(nat),atoms%jri(nat),dpdot_re)
         CALL intgr3(dpj_im,atoms%rmsh(1,nat),atoms%dx(nat),atoms%jri(nat),dpdot_im)
         tmt = tmt + CMPLX(dpdot_re,dpdot_im)*atoms%neq(nat)
      END DO

      mt_int = mt_int + tmt

   END SUBROUTINE dfpt_int_mt

   SUBROUTINE dfpt_int_mt_sf(atoms, sphhar, sym, iDir, nat, mt_conj, mt_pure, mt_pure_im, sf_int)
      USE m_gaunt, ONLY: gaunt1

      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_sphhar), INTENT(IN) :: sphhar
      TYPE(t_sym),    INTENT(IN) :: sym

      INTEGER, INTENT(IN) :: iDir, nat

      REAL, INTENT(IN)  :: mt_conj(:,0:,:), mt_pure(:,0:,:), mt_pure_im(:,0:,:)

      COMPLEX, INTENT(INOUT) :: sf_int

      INTEGER :: lhPr, lh, lPr, l, iMemPr, iMem, mPr, m, mPrPr, nd
      REAL    :: f_prod_re, f_prod_im, gaunt_coef
      COMPLEX :: f_prod, tmt

      COMPLEX :: trafo_mat(3,-1:1)

      nd = sym%ntypsy(nat)

      trafo_mat = 0.0

      trafo_mat(1,-1) =  1
      trafo_mat(1, 1) = -1
      trafo_mat(2,-1) = ImagUnit
      trafo_mat(2, 1) = ImagUnit
      trafo_mat(3, 0) = SQRT(2.0)

      trafo_mat = SQRT(2*pi_const/3) * trafo_mat

      tmt = CMPLX(0.0,0.0)
      DO lhPr = 0, sphhar%nlh(sym%ntypsy(nat))
         lPr = sphhar%llh(lhPr,nd)
         DO lh = 0, sphhar%nlh(sym%ntypsy(nat))
            l = sphhar%llh(lh,nd)
            f_prod_re =  mt_pure(atoms%jri(nat),lh,nat)    * mt_conj(atoms%jri(nat),lhPr,nat)
            f_prod_im =  mt_pure_im(atoms%jri(nat),lh,nat) * mt_conj(atoms%jri(nat),lhPr,nat)
            f_prod = CMPLX(f_prod_re,f_prod_im)
            DO iMemPr = 1, sphhar%nmem(lhPr,nd)
               mPr = sphhar%mlh(iMemPr,lhPr,nd)
               DO iMem = 1, sphhar%nmem(lh,nd)
                  m = sphhar%mlh(iMem,lh,nd)
                  IF (lPr<ABS(l-1).OR.lPr>(l+1)) CYCLE ! |l-1| <= l' <= l+1
                  IF (MOD(lPr+l,2)/=1) CYCLE ! l' + l + 1 even
                  mPrPr = mPr - m ! m' = m'' + m
                  IF (ABS(mPrPr)>1) CYCLE ! |m''| <= 1
                  gaunt_coef = gaunt1(lPr,1,l,mPr,mPrPr,m,atoms%lmaxd)
                  tmt = tmt + f_prod * CONJG(sphhar%clnu(iMemPr,lhPr,nd)) * sphhar%clnu(iMem,lh,nd) &
                            * gaunt_coef * trafo_mat(iDir,mPrPr)
               END DO
            END DO
         END DO
      END DO

      sf_int = sf_int + tmt

   END SUBROUTINE dfpt_int_mt_sf

   SUBROUTINE dfpt_dynmat_eigen(fi, results, results1, xcpot, fmpi, mpdata, hybdat, enpara, nococonv, &
                                stars, starsq, sphhar, inden, hub1data, vx, v, v1real, v1imag, &
                                eig_id, dfpt_eig_id, dfpt_eig_id2, iDir_col, iDtype_col, iDir_row, iDtype_row, &
                                theta1_pw0, theta1_pw, bqpt, l_real, eigen_term, killcont, q_eig_id)

      USE m_types
      USE m_constants
      USE m_eigen_hssetup
      USE m_pot_io
      USE m_eigen_diag
      USE m_local_hamiltonian
      USE m_util
      USE m_eig66_io, ONLY : write_eig, read_eig
      USE m_xmlOutput
      USE m_types_mpimat
      USE m_dfpt_tlmplm
      USE m_npy

! TODO: One brigth day, these things will also be relevant for DFPT.
!       We cannot keep doing small systems on small CPUs forever.
!#ifdef _OPENACC
!         USE cublas
!#define CPP_zgemv cublaszgemv
!#else
#define CPP_zgemv zgemv
!#endif

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_results),INTENT(INOUT):: results, results1
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_mpdata), intent(inout):: mpdata
      TYPE(t_hybdat), INTENT(INOUT):: hybdat
      TYPE(t_enpara),INTENT(INOUT) :: enpara
      TYPE(t_nococonv),INTENT(IN)  :: nococonv
      TYPE(t_stars),INTENT(IN)     :: stars, starsq
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_potden),INTENT(IN)    :: inden !
      TYPE(t_hub1data),INTENT(INOUT):: hub1data
      TYPE(t_potden), INTENT(IN)   :: vx
      TYPE(t_potden),INTENT(IN)    :: v, v1real, v1imag

      ! Scalar Arguments
      INTEGER, INTENT(IN)    :: eig_id, dfpt_eig_id, iDir_col, iDtype_col, iDir_row, iDtype_row, dfpt_eig_id2
      COMPLEX,            INTENT(IN)     :: theta1_pw0(:), theta1_pw(:)

      REAL,    INTENT(IN) :: bqpt(3)
      LOGICAL, INTENT(IN) :: l_real

      COMPLEX, INTENT(INOUT) :: eigen_term

      INTEGER, OPTIONAL, INTENT(IN) :: killcont(11), q_eig_id

      ! Local Scalars
      INTEGER jsp,nk,ne_all,ne_found,neigd2
      INTEGER nk_i,n_size,n_rank
      INTEGER err
      REAL :: q_loop(3)
      ! Local Arrays
      INTEGER              :: ierr, nbands, nbands1, nbasfcn, nbasfcnq, noccbd, nbandsq, nbands2 

      REAL,    ALLOCATABLE :: bkpt(:)
      REAL,    ALLOCATABLE :: eig(:), eig1(:), we(:), we1(:)

      TYPE(t_tlmplm)            :: tdV1, tdmod, td
      TYPE(t_usdus)             :: ud, uddummy
      TYPE(t_lapw)              :: lapw, lapwq
      TYPE(t_kpts)              :: kpts_mod !kqpts ! basically kpts, but with q added onto each one.
      TYPE(t_hub1data)          :: hub1datadummy
      TYPE (t_mat)              :: zMat, zMat1, zMatq, zMat2
      CLASS(t_mat), ALLOCATABLE :: hmat1,smat1,hmat1q,smat1q,hmat2,smat2,vmat2

      ! Variables for HF or fi%hybinp functional calculation
      INTEGER                   :: comm(fi%kpts%nkpt),irank2(fi%kpts%nkpt),isize2(fi%kpts%nkpt), dealloc_stat
      character(len=300)        :: errmsg

      INTEGER :: iEig, ikGq, iqdir
      INTEGER :: imlo, ilo, iklo, l, ikg, ikglo !!! Test!
      REAL :: gext(3)
      COMPLEX :: we_loop, we1_loop, eig_loop, eig1_loop, eigen_s1q, eigen_h1qs1q, eigen_bonus
      COMPLEX :: eigen_e1, eigen_w1e1, eigen_s2w1e1, eigen_h2s2w1e1, eigen_w1h2s2w1e1

      COMPLEX, ALLOCATABLE :: tempVec(:), tempVecq(:), z_loop(:), z1_loop(:), ztest_loop(:), zq_loop(:)

      INTEGER, ALLOCATABLE :: k_selection(:)

      COMPLEX  zdotc
      EXTERNAL zdotc
      
      eigen_e1 = CMPLX(0.0,0.0)
      eigen_w1e1 = CMPLX(0.0,0.0)
      eigen_s2w1e1 = CMPLX(0.0,0.0)
      eigen_h2s2w1e1 = CMPLX(0.0,0.0)
      eigen_w1h2s2w1e1 = CMPLX(0.0,0.0)
      eigen_s1q = CMPLX(0.0,0.0)
      eigen_h1qs1q = CMPLX(0.0,0.0)
      eigen_bonus = CMPLX(0.0,0.0) 
      ALLOCATE(k_selection(16))
      k_selection = [25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
      !kqpts = fi%kpts

      kpts_mod = fi%kpts
      ! Modify this from kpts only in DFPT case.
      ALLOCATE(bkpt(3))
      DO nk_i = 1, size(fmpi%k_list)
         !kqpts%bk(:, nk_i) = kqpts%bk(:, nk_i) + bqpt
         nk=fmpi%k_list(nk_i)
         bkpt = fi%kpts%bk(:, nk)
         DO iqdir = 1, 3
            !IF (bkpt(iqdir)+bqpt(iqdir)>=0.5) bkpt(iqdir) = bkpt(iqdir) - 1.0
            !IF (bkpt(iqdir)+bqpt(iqdir)<-0.5) bkpt(iqdir) = bkpt(iqdir) + 1.0
            !IF (bkpt(iqdir)+bqpt(iqdir)>=0.5.AND.ABS(bqpt(iqdir))>1e-8) bkpt(iqdir) = bkpt(iqdir) - 1.0
            !IF (bkpt(iqdir)+bqpt(iqdir)<-0.5.AND.ABS(bqpt(iqdir))>1e-8) bkpt(iqdir) = bkpt(iqdir) + 1.0
         END DO
         kpts_mod%bk(:, nk) = bkpt
      END DO

      call ud%init(fi%atoms,fi%input%jspins)
      call uddummy%init(fi%atoms,fi%input%jspins)
      ALLOCATE(eig(fi%input%neig))
    
      CALL local_ham(sphhar,fi%atoms,fi%sym,fi%noco,nococonv,enpara,fmpi,v,vx,inden,fi%input,fi%hub1inp,hub1data,td,ud,0.0)
      ! Get matrix elements of perturbed potential and modified H/S in DFPT case.
      hub1datadummy = hub1data

      CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,v,fmpi,tdV1,v1real,v1imag,.FALSE.,iDtype_col)

      CALL local_ham(sphhar,fi%atoms,fi%sym,fi%noco,nococonv,enpara,fmpi,v,vx,inden,fi%input,fi%hub1inp,hub1datadummy,tdmod,uddummy,0.0,.true.)

      DO jsp = MERGE(1,1,fi%noco%l_noco), MERGE(1,fi%input%jspins,fi%noco%l_noco)
         k_loop:DO nk_i = 1,size(fmpi%k_list)
            nk = fmpi%k_list(nk_i)
            bkpt = kpts_mod%bk(:, nk)
            q_loop = bqpt

            CALL lapw%init(fi%input,fi%noco,nococonv,fi%kpts,fi%atoms,fi%sym,nk,fi%cell,fmpi)
            CALL lapwq%init(fi%input,fi%noco,nococonv,kpts_mod,fi%atoms,fi%sym,nk,fi%cell,fmpi,q_loop)

            we  = results%w_iks(:,nk,jsp)
            we1 = results1%w_iks(:,nk,jsp)
            eig = results%eig(:,nk,jsp)
            eig1 = results1%eig(:,nk,jsp)

            noccbd = COUNT(we>1.e-8)

            nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
            nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

            CALL zMat%init(l_real,nbasfcn,noccbd)
            CALL zMat1%init(.FALSE.,nbasfcnq,noccbd)
            CALL zMat2%init(.FALSE.,nbasfcnq,noccbd)
            IF (PRESENT(q_eig_id)) CALL zMatq%init(l_real,nbasfcnq,noccbd)

            ALLOCATE(tempVec(nbasfcn))
            IF (PRESENT(q_eig_id)) ALLOCATE(tempVecq(nbasfcnq))
            ALLOCATE(z_loop(nbasfcn),z1_loop(nbasfcnq))
            ALLOCATE(ztest_loop(nbasfcn))
            IF (PRESENT(q_eig_id)) ALLOCATE(zq_loop(nbasfcnq))

            CALL read_eig(eig_id,      nk,jsp,neig=nbands,zmat=zMat)
            CALL read_eig(dfpt_eig_id, nk,jsp,neig=nbands1,zmat=zMat1)
            CALL read_eig(dfpt_eig_id2,nk,jsp,neig=nbands2,zmat=zMat2)
            IF (PRESENT(q_eig_id)) CALL read_eig(q_eig_id, nk,jsp,neig=nbandsq,zmat=zMatq)

            !!!! Test:
            !DO ikG = 1, lapw%nv(jsp)
            !   gExt = MATMUL(fi%cell%bmat,lapw%vk(:, ikG, jsp))
            !   IF (zMat%l_real) THEN
            !      zMat1%data_c(ikG,:) = -ImagUnit * gExt(idir_row) * zMat%data_r(ikG, :noccbd)
            !   ELSE
            !      zMat1%data_c(ikG,:) = -ImagUnit * gExt(idir_row) * zMat%data_c(ikG, :noccbd)
            !   END IF
            !END DO

            !DO ikG = lapw%nv(jsp) + 1, lapw%nv(jsp) + fi%atoms%nlo(iDtype_row)
            !   iLo = ikG-lapw%nv(jsp)
            !   l = fi%atoms%llo(iLo, iDtype_row)
            !   DO imLo = 1, 2*l+1
            !      ikLo = lapw%kvec(imLo,iLo,iDtype_row)
            !      ikGLo = lapw%nv(jsp) + lapw%index_lo(iLo,iDtype_row) + imLo
            !         gExt = MATMUL(fi%cell%bmat,lapw%vk(:,ikLo, jsp))
            !         IF (zMat%l_real) THEN
            !            zMat1%data_c(ikGLo,:) = -ImagUnit * gExt(idir_row) * zMat%data_r(ikGLo, :)
            !         ELSE
            !            zMat1%data_c(ikGLo,:) = -ImagUnit * gExt(idir_row) * zMat%data_c(ikGLo, :)
            !         END IF
            !   END DO
            !   gExt = MATMUL(fi%cell%bmat,lapw%vk(:, ikG, jsp))
            !END DO

            CALL timestart("Setup of H&S matrices")
            IF (.NOT.PRESENT(q_eig_id)) THEN
               CALL dfpt_dynmat_hssetup(jsp, fmpi, fi, enpara, nococonv, starsq, stars, &
                                        ud, tdmod, tdV1, lapw, lapwq, iDir_row, iDtype_row, iDir_col, iDtype_col, theta1_pw0, theta1_pw, &
                                        smat1, hmat1, smat1q, hmat1q, smat2, hmat2, nk, killcont)
            ELSE
               CALL dfpt_dynmat_hssetup(jsp, fmpi, fi, enpara, nococonv, starsq, stars, &
                                        ud, tdmod, tdV1, lapw, lapwq, iDir_row, iDtype_row, iDir_col, iDtype_col, theta1_pw0, theta1_pw, &
                                        smat1, hmat1, smat1q, hmat1q, smat2, hmat2, nk, killcont, vmat2)
            END IF
            CALL timestop("Setup of H&S matrices")

            DO iEig = 1, noccbd
               eig_loop  = eig(iEig)
               eig1_loop = eig1(iEig)
               we_loop   = (2.0/fi%input%jspins)*we(iEig)
               we1_loop  = (2.0/fi%input%jspins)*we1(iEig)
               IF (l_real) THEN
                  z_loop = CMPLX(1.0,0.0)*zMat%data_r(:,iEig)
                  IF (PRESENT(q_eig_id)) zq_loop = CMPLX(1.0,0.0)*zMatq%data_r(:,iEig)
               ELSE
                  z_loop    = zMat%data_c(:,iEig)
                  IF (PRESENT(q_eig_id)) zq_loop = zMatq%data_c(:,iEig)
               END IF

               z1_loop = zMat1%data_c(:,iEig)

               CALL CPP_zgemv('N',nbasfcn,nbasfcn,-we_loop*eig1_loop,smat1%data_c,nbasfcn,z_loop,1,CMPLX(0.0,0.0),tempVec,1)
               eigen_e1 = eigen_e1 + zdotc(nbasfcn,z_loop,1,tempVec,1)
               CALL CPP_zgemv('N',nbasfcn,nbasfcn,-we1_loop*eig_loop,smat1%data_c,nbasfcn,z_loop,1,CMPLX(1.0,0.0),tempVec,1)
               eigen_w1e1 = eigen_w1e1 + zdotc(nbasfcn,z_loop,1,tempVec,1)
               CALL CPP_zgemv('N',nbasfcn,nbasfcn,-we_loop*eig_loop,smat2%data_c,nbasfcn,z_loop,1,CMPLX(1.0,0.0),tempVec,1)
               eigen_s2w1e1 = eigen_s2w1e1 + zdotc(nbasfcn,z_loop,1,tempVec,1)
               CALL CPP_zgemv('N',nbasfcn,nbasfcn,we_loop,hmat2%data_c,nbasfcn,z_loop,1,CMPLX(1.0,0.0),tempVec,1)
               eigen_h2s2w1e1 = eigen_h2s2w1e1 + zdotc(nbasfcn,z_loop,1,tempVec,1)
               CALL CPP_zgemv('N',nbasfcn,nbasfcn,we1_loop,hmat1%data_c,nbasfcn,z_loop,1,CMPLX(1.0,0.0),tempVec,1)
               eigen_w1h2s2w1e1 = eigen_w1h2s2w1e1 + zdotc(nbasfcn,z_loop,1,tempVec,1)

               eigen_term = eigen_term + zdotc(nbasfcn,z_loop,1,tempVec,1)

               IF (PRESENT(q_eig_id)) THEN
                  CALL CPP_zgemv('N',nbasfcnq,nbasfcn,we_loop,vmat2%data_c,nbasfcnq,z_loop,1,CMPLX(0.0,0.0),tempVecq,1)
                  eigen_term = eigen_term + zdotc(nbasfcnq,zq_loop,1,tempVecq,1)
               END IF

               CALL CPP_zgemv('C',nbasfcnq,nbasfcn,-we_loop*eig_loop,smat1q%data_c,nbasfcnq,z1_loop,1,CMPLX(0.0,0.0),tempVec,1)
               eigen_s1q = eigen_s1q + 2.0*zdotc(nbasfcn,z_loop,1,tempVec,1)
               CALL CPP_zgemv('C',nbasfcnq,nbasfcn,we_loop,hmat1q%data_c,nbasfcnq,z1_loop,1,CMPLX(1.0,0.0),tempVec,1)
               eigen_h1qs1q = eigen_h1qs1q + 2.0*zdotc(nbasfcn,z_loop,1,tempVec,1)

               eigen_term = eigen_term + 2.0*zdotc(nbasfcn,z_loop,1,tempVec,1)
               
               ! Correction:
               CALL CPP_zgemv('C',nbasfcnq,nbasfcn,we_loop,smat1q%data_c,nbasfcnq,zMat2%data_c(:,iEig),1,CMPLX(0.0,0.0),tempVec,1)
               eigen_bonus = eigen_bonus + 2.0*zdotc(nbasfcn,z_loop,1,tempVec,1)
               eigen_term = eigen_term + 2.0 * zdotc(nbasfcn,z_loop,1,tempVec,1)
            END DO

            DEALLOCATE(tempVec)
            IF (PRESENT(q_eig_id)) DEALLOCATE(tempVecq)
            DEALLOCATE(z_loop,z1_loop)
            IF (PRESENT(q_eig_id)) DEALLOCATE(zq_loop) 
            DEALLOCATE(ztest_loop)
            CALL smat1%free()
            CALL hmat1%free()
            CALL smat1q%free()
            CALL hmat1q%free()
            CALL smat2%free()
            CALL hmat2%free()
            IF (PRESENT(q_eig_id)) CALL vmat2%free()
            DEALLOCATE(hmat1,smat1,hmat1q,smat1q,hmat2,smat2, stat=dealloc_stat, errmsg=errmsg)
            IF (PRESENT(q_eig_id)) DEALLOCATE(vmat2)
            if(dealloc_stat /= 0) call juDFT_error("deallocate failed one of the matrices",&
                                                   hint=errmsg, calledby="dfpt_dynmat.F90")

            ! Output results
            CALL timestart("EV output")

#if defined(CPP_MPI)
            ! RMA synchronization
            CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif
            CALL timestop("EV output")

            !IF (allocated(zmat)) THEN
             CALL zMat%free()
             CALL zMat1%free()
             CALL zMat2%free()
              !deallocate(zMat)
            !ENDIF
         END DO  k_loop
      END DO ! spin loop ends
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_term,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_e1,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_w1e1,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_s2w1e1,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_h2s2w1e1,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_w1h2s2w1e1,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_s1q,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_h1qs1q,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,eigen_bonus,1,MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
#endif
      IF (fmpi%irank==0) THEN
         write(9989,*) "eigen -w0e1s0 :", eigen_e1
         write(9989,*) "eigen -w1e0s0 :", eigen_w1e1-eigen_e1
         write(9989,*) "eigen -w0e0s2 :", eigen_s2w1e1-eigen_w1e1-eigen_e1
         write(9989,*) "eigen w0h2    :", eigen_h2s2w1e1-eigen_s2w1e1-eigen_w1e1-eigen_e1
         write(9989,*) "eigen w1h1    :", eigen_w1h2s2w1e1-eigen_h2s2w1e1-eigen_s2w1e1-eigen_w1e1-eigen_e1
         write(9989,*) "eigen -w0e0s1q:", eigen_s1q
         write(9989,*) "eigen w0h1q   :", eigen_h1qs1q - eigen_s1q
         write(9989,*) "eigen corr    :", eigen_bonus
      END IF
   END SUBROUTINE

   SUBROUTINE dfpt_dynmat_hssetup(isp, fmpi, fi, enpara, nococonv, starsq, stars, &
                            ud, td, tdV1, lapw, lapwq, iDir_row, iDtype_row, iDir_col, iDtype_col, theta1_pw0, theta1_pw, &
                            smat1_final, hmat1_final, smat1q_final, hmat1q_final, smat2_final, hmat2_final, nk, killcont, vmat2_final)
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
      TYPE(t_stars),      INTENT(IN)     :: starsq, stars
      TYPE(t_enpara),     INTENT(IN)     :: enpara
      TYPE(t_nococonv),   INTENT(IN)     :: nococonv
      TYPE(t_usdus),      INTENT(IN)     :: ud
      TYPE(t_tlmplm),     INTENT(IN)     :: td, tdV1
      TYPE(t_lapw),       INTENT(IN)     :: lapw, lapwq
      INTEGER,            INTENT(IN)     :: iDir_row, iDtype_row, iDir_col, iDtype_col
      COMPLEX,            INTENT(IN)     :: theta1_pw0(:), theta1_pw(:)
      CLASS(t_mat), ALLOCATABLE, INTENT(INOUT)   :: smat1_final, hmat1_final, smat1q_final, hmat1q_final, smat2_final, hmat2_final
      CLASS(t_mat), OPTIONAL, ALLOCATABLE, INTENT(INOUT)   :: vmat2_final
      INTEGER,      INTENT(IN)     :: nk, killcont(11)

      CLASS(t_mat), ALLOCATABLE :: smat1(:, :), hmat1(:, :), smat1q(:, :), hmat1q(:, :), smat2(:, :), hmat2(:, :), vmat2(:,:)

      INTEGER :: i, j, nspins

      nspins = MERGE(2, 1, fi%noco%l_noco)
      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::smat1(nspins, nspins), hmat1(nspins, nspins))
         ALLOCATE (t_mat::smat1q(nspins, nspins), hmat1q(nspins, nspins))
         ALLOCATE (t_mat::smat2(nspins, nspins), hmat2(nspins, nspins))
         IF (PRESENT(vmat2_final)) THEN
            ALLOCATE (t_mat::vmat2(nspins, nspins))
         END IF
      ELSE
         ALLOCATE (t_mpimat::smat1(nspins, nspins), hmat1(nspins, nspins))
         ALLOCATE (t_mpimat::smat1q(nspins, nspins), hmat1q(nspins, nspins))
         ALLOCATE (t_mpimat::smat2(nspins, nspins), hmat2(nspins, nspins))
         IF (PRESENT(vmat2_final)) THEN
            ALLOCATE (t_mpimat::vmat2(nspins, nspins))
         END IF 
      END IF

      DO i = 1, nspins
         DO j = 1, nspins
            CALL smat1(i, j)%init(.FALSE., lapw%nv(i) + fi%atoms%nlotot, lapw%nv(j) + fi%atoms%nlotot, fmpi%sub_comm, .false.)
            CALL hmat1(i, j)%init(smat1(i, j))
            CALL smat1q(i, j)%init(.FALSE., lapwq%nv(i) + fi%atoms%nlotot, lapw%nv(j) + fi%atoms%nlotot, fmpi%sub_comm, .false.)
            CALL hmat1q(i, j)%init(smat1q(i, j))
            CALL smat2(i, j)%init(.FALSE., lapw%nv(i) + fi%atoms%nlotot, lapw%nv(j) + fi%atoms%nlotot, fmpi%sub_comm, .false.)
            CALL hmat2(i, j)%init(smat2(i, j))
            IF (PRESENT(vmat2_final)) THEN
               CALL vmat2(i, j)%init(.FALSE., lapwq%nv(i) + fi%atoms%nlotot, lapw%nv(j) + fi%atoms%nlotot, fmpi%sub_comm, .false.)
            END IF
         END DO
      END DO

      CALL timestart("Interstitial part")
      CALL dfpt_dynmat_hs_int(fi%noco, starsq, stars, lapwq, lapw, fmpi, fi%cell%bbmat, isp, theta1_pw0, theta1_pw, &
                              smat1, hmat1, smat1q, hmat1q, killcont(1:4))
      CALL timestop("Interstitial part")

      CALL timestart("MT part")
      DO i = 1, nspins; DO j = 1, nspins
            !!$acc enter data copyin(hmat(i,j),smat(i,j))
            !!$acc enter data copyin(hmat(i,j)%data_r,smat(i,j)%data_r,hmat(i,j)%data_c,smat(i,j)%data_c)
      END DO; END DO
      IF (.NOT.PRESENT(vmat2_final)) THEN
         CALL dfpt_dynmat_hsmt(fi%atoms, fi%sym, enpara, isp, iDir_row, iDtype_row, iDir_col, iDtype_col, fi%input, fmpi, fi%noco, nococonv, fi%cell, &
                               lapw, lapwq, ud, td, tdV1, hmat1, smat1, hmat1q, smat1q, hmat2, smat2, nk, killcont(5:11))
      ELSE
         CALL dfpt_dynmat_hsmt(fi%atoms, fi%sym, enpara, isp, iDir_row, iDtype_row, iDir_col, iDtype_col, fi%input, fmpi, fi%noco, nococonv, fi%cell, &
                               lapw, lapwq, ud, td, tdV1, hmat1, smat1, hmat1q, smat1q, hmat2, smat2, nk, killcont(5:11), vmat2)
      END IF
      DO i = 1, nspins; DO j = 1, nspins; if (hmat1(1, 1)%l_real) THEN
            !!$acc exit data copyout(hmat(i,j)%data_r,smat(i,j)%data_r) delete(hmat(i,j)%data_c,smat(i,j)%data_c)
            !!$acc exist data delete(hmat(i,j),smat(i,j))
         ELSE
            !!$acc exit data copyout(hmat(i,j)%data_c,smat(i,j)%data_c) delete(hmat(i,j)%data_r,smat(i,j)%data_r)
            !!$acc exit data delete(hmat(i,j),smat(i,j))
         END IF; END DO; END DO
      CALL timestop("MT part")

      !Now copy the data into final matrix
      ! Collect the four fi%noco parts into a single matrix
      ! In collinear case only a copy is done
      ! In the parallel case also a redistribution happens
      ALLOCATE (smat1_final, mold=smat1(1, 1))
      ALLOCATE (hmat1_final, mold=smat1(1, 1))
      ALLOCATE (smat1q_final, mold=smat1q(1, 1))
      ALLOCATE (hmat1q_final, mold=smat1q(1, 1))
      ALLOCATE (smat2_final, mold=smat2(1, 1))
      ALLOCATE (hmat2_final, mold=smat2(1, 1))
      IF (PRESENT(vmat2_final)) ALLOCATE (vmat2_final, mold=vmat2(1, 1))

      CALL timestart("Matrix redistribution")
      CALL dfpt_eigen_redist_matrix(fmpi, lapw, lapw, fi%atoms, smat1, smat1_final)
      CALL dfpt_eigen_redist_matrix(fmpi, lapw, lapw, fi%atoms, hmat1, hmat1_final, smat1_final)
      CALL dfpt_eigen_redist_matrix(fmpi, lapwq, lapw, fi%atoms, smat1q, smat1q_final)
      CALL dfpt_eigen_redist_matrix(fmpi, lapwq, lapw, fi%atoms, hmat1q, hmat1q_final, smat1q_final)
      CALL dfpt_eigen_redist_matrix(fmpi, lapw, lapw, fi%atoms, smat2, smat2_final)
      CALL dfpt_eigen_redist_matrix(fmpi, lapw, lapw, fi%atoms, hmat2, hmat2_final, smat2_final)
      IF (PRESENT(vmat2_final)) CALL dfpt_eigen_redist_matrix(fmpi, lapw, lapw, fi%atoms, vmat2, vmat2_final)
      CALL timestop("Matrix redistribution")
   END SUBROUTINE
END MODULE m_dfpt_dynmat
