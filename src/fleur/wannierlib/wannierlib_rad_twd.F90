!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_rad_twd
  USE m_juDFT
  USE m_constants
  USE m_intgr, ONLY: intgr3
  USE m_types_atoms
  USE m_types_usdus
  USE m_types_radfun
  USE m_clebsch
  USE m_types_wannierlib
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_rad_twd(wannierlib, atoms, nwfs, ikpt, usdus, radfun,  jspin, rads)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: wannierlib
    TYPE(t_atoms), INTENT(IN) :: atoms
    INTEGER, INTENT(IN) :: nwfs
    INTEGER, INTENT(IN) :: ikpt
    TYPE(t_usdus), INTENT(IN) :: usdus
    TYPE(t_radfun), INTENT(IN) :: radfun(atoms%ntype)
    INTEGER, INTENT(IN) :: jspin
    REAL, INTENT(OUT) :: rads(nwfs, 0:3, atoms%jmtd, 2)

    INTEGER :: nwf, ntyp, j, l, lo
    REAL :: alpha, wronk, radi
    REAL :: acft, bcft, aa, bb, a1, b1, rho
    REAL :: radf(atoms%jmtd)

    

    CALL timestart('wannierlib_rad_twd')

    DO nwf = 1, nwfs
      IF (ABS(wannierlib%proj_regio(nwf) - 1.0) >= 1.0e-5) THEN
        CALL juDFT_error('no projections outside muffin tins', calledby='wannierlib_rad_twd')
      END IF

      ntyp = wannierlib%proj_ntype(nwf)
      alpha = wannierlib%proj_zona(nwf)
      rads(nwf, :, :, :) = 0.0

      IF (wannierlib%proj_rwf(nwf) == -5) THEN
        aa = 2.0 * (alpha ** 1.5)
        bb = alpha
        a1 = aa * EXP(-bb * atoms%rmsh(atoms%jri(ntyp), ntyp))
        b1 = -aa * bb * EXP(-bb * atoms%rmsh(atoms%jri(ntyp), ntyp))

        DO l = 0, 3
          wronk = usdus%us(l, ntyp, jspin) * usdus%duds(l, ntyp, jspin) - &
                  usdus%uds(l, ntyp, jspin) * usdus%dus(l, ntyp, jspin)
          acft = (a1 * usdus%duds(l, ntyp, jspin) - b1 * usdus%uds(l, ntyp, jspin)) / wronk
          bcft = (b1 * usdus%us(l, ntyp, jspin) - a1 * usdus%dus(l, ntyp, jspin)) / wronk
          DO j = 1, atoms%jri(ntyp)
            rads(nwf, l, j, :) = acft * radfun(ntyp)%r(j, :, 1, l, jspin) + &
                                 bcft * radfun(ntyp)%r(j, :, 2, l, jspin)
          END DO
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == 0) THEN
        DO l = 0, 3
          DO j = 1, atoms%jri(ntyp)
            rads(nwf, l, j, :) = radfun(ntyp)%r(j, :, 1, l, jspin) * (1.0 - ABS(alpha)) + &
                                 alpha * radfun(ntyp)%r(j, :, 2, l, jspin)
          END DO
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == -6) THEN
        IF (atoms%nlod < 1) CALL juDFT_error('nlod<1', calledby='wannierlib_rad_twd')
        DO l = 0, 3
          DO j = 1, atoms%jri(ntyp)
            rads(nwf, l, j, :) = radfun(ntyp)%r(j, :, 3, l, jspin)
          END DO
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == -7) THEN
        IF (atoms%nlod < 2) CALL juDFT_error('nlod<2', calledby='wannierlib_rad_twd')
        DO l = 0, 3
          DO j = 1, atoms%jri(ntyp)
            rads(nwf, l, j, :) = radfun(ntyp)%r(j, :, 4, l, jspin)
          END DO
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == -8) THEN
        IF (atoms%nlod < 3) CALL juDFT_error('nlod<3', calledby='wannierlib_rad_twd')
        DO l = 0, 3
          DO j = 1, atoms%jri(ntyp)
            rads(nwf, l, j, :) = radfun(ntyp)%r(j, :, 5, l, jspin)
          END DO
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == -9) THEN
        rads(nwf, :, :, :) = 0.0
        DO lo = 1, atoms%nlo(ntyp)
          l = atoms%llo(lo, ntyp)
          DO j = 1, atoms%jri(ntyp)
            rads(nwf, l, j, :) = radfun(ntyp)%r(j, :, 6, l, jspin)
          END DO
        END DO
      ELSE IF ((wannierlib%proj_rwf(nwf) > -5) .AND. (wannierlib%proj_rwf(nwf) < 0)) THEN
        DO l = 0, 3
          DO j = 1, atoms%jri(ntyp)
            rads(nwf, l, j, :) = radfun(ntyp)%r(j, :, 1, ABS(wannierlib%proj_rwf(nwf)) - 1, jspin) * (1.0 - ABS(alpha)) + &
                                 alpha * radfun(ntyp)%r(j, :, 2, ABS(wannierlib%proj_rwf(nwf)) - 1, jspin)
          END DO
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == 1) THEN
        DO j = 1, atoms%jri(ntyp)
          rho = 2.0 * alpha * atoms%rmsh(j, ntyp)
          rads(nwf, 0, j, 1) = 2.0 * (alpha ** 1.5) * EXP(-rho / 2.0)
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == 2) THEN
        DO j = 1, atoms%jri(ntyp)
          rho = alpha * atoms%rmsh(j, ntyp)
          rads(nwf, 0, j, 1) = (1.0 / (2.0 * SQRT(2.0))) * (2.0 - rho) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 1, j, 1) = (1.0 / (2.0 * SQRT(6.0))) * rho * (alpha ** 1.5) * EXP(-rho / 2.0)
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == 3) THEN
        DO j = 1, atoms%jri(ntyp)
          rho = 2.0 * alpha * atoms%rmsh(j, ntyp) / 3.0
          rads(nwf, 0, j, 1) = (1.0 / (9.0 * SQRT(3.0))) * (6.0 - 6.0 * rho + rho ** 2) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 1, j, 1) = (1.0 / (9.0 * SQRT(6.0))) * (4.0 - rho) * rho * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 2, j, 1) = (1.0 / (9.0 * SQRT(30.0))) * rho * rho * (alpha ** 1.5) * EXP(-rho / 2.0)
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == 4) THEN
        DO j = 1, atoms%jri(ntyp)
          rho = 2.0 * alpha * atoms%rmsh(j, ntyp) / 4.0
          rads(nwf, 0, j, 1) = (1.0 / 96.0) * (24.0 - 36.0 * rho + 12.0 * rho ** 2 - rho ** 3) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 1, j, 1) = (1.0 / (32.0 * SQRT(15.0))) * rho * (20.0 - 10.0 * rho + rho ** 2) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 2, j, 1) = (1.0 / (96.0 * SQRT(5.0))) * rho * rho * (6.0 - rho) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 3, j, 1) = (1.0 / (96.0 * SQRT(35.0))) * rho * rho * rho * (alpha ** 1.5) * EXP(-rho / 2.0)
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == 5) THEN
        DO j = 1, atoms%jri(ntyp)
          rho = 2.0 * alpha * atoms%rmsh(j, ntyp) / 5.0
          rads(nwf, 0, j, 1) = (1.0 / (300.0 * SQRT(5.0))) * (120.0 - 240.0 * rho + 120.0 * rho ** 2 - 20.0 * rho ** 3 + rho ** 4) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 1, j, 1) = (1.0 / (150.0 * SQRT(30.0))) * rho * (120.0 - 90.0 * rho + 18.0 * rho ** 2 - rho ** 3) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 2, j, 1) = (1.0 / (150.0 * SQRT(70.0))) * rho * rho * (42.0 - 14.0 * rho + rho * rho) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 3, j, 1) = (1.0 / (300.0 * SQRT(70.0))) * rho * rho * rho * (8.0 - rho) * (alpha ** 1.5) * EXP(-rho / 2.0)
        END DO
      ELSE IF (wannierlib%proj_rwf(nwf) == 6) THEN
        DO j = 1, atoms%jri(ntyp)
          rho = 2.0 * alpha * atoms%rmsh(j, ntyp) / 6.0
          rads(nwf, 0, j, 1) = (1.0 / 2160.0 * SQRT(6.0)) * (720.0 - 1800.0 * rho + 1200.0 * rho * rho - 300.0 * rho ** 3 + 30.0 * rho ** 4 - rho ** 5) * (alpha ** 1.5) * EXP(-rho / 2.0)
          rads(nwf, 1, j, 1) = (1.0 / 432.0 * SQRT(210.0)) * rho * (840.0 - 840.0 * rho + 252.0 * rho * rho - 28.0 * rho ** 3 + rho ** 4) * (alpha ** 1.5) * EXP(-rho / 2.0)
        END DO
      ELSE
        CALL juDFT_error('radial function is not tabulated', calledby='wannierlib_rad_twd')
      END IF

      IF (ikpt == 1) THEN
        DO l = 0, 3
          DO j = 1, atoms%jri(ntyp)
            radf(j) = atoms%rmsh(j, ntyp) * atoms%rmsh(j, ntyp) * rads(nwf, l, j, 1) * rads(nwf, l, j, 1)
          END DO
          CALL intgr3(radf, atoms%rmsh(1, ntyp), atoms%dx(ntyp), atoms%jri(ntyp), radi)
          WRITE (oUnit, *)
          WRITE (oUnit, *) 'Wannier Function N:', nwf
          WRITE (oUnit, *) 'angular momentum', l
          WRITE (oUnit, *) 'radial function at the MT boundary:', rads(nwf, l, atoms%jri(ntyp), 1)
          WRITE (oUnit, *) 'norma =', radi
        END DO
      END IF
    END DO

    CALL timestop('wannierlib_rad_twd')
  END SUBROUTINE wannierlib_rad_twd

 SUBROUTINE wannierlib_soc_tlmw(nwfs, lwf, jwf, jmwf, jspin, tlmwf)
    INTEGER, INTENT(IN) :: nwfs
    INTEGER, INTENT(IN) :: lwf(:)
    REAL, INTENT(IN) :: jwf(:), jmwf(:)
    INTEGER, INTENT(IN) :: jspin
    COMPLEX, INTENT(OUT) :: tlmwf(0:3, -3:3, nwfs)

    INTEGER :: nwf, l, m
    REAL :: j, jm

    CALL timestart('wannierlib_soc_tlmw')
    DO nwf = 1, nwfs
      l = lwf(nwf)
      IF (l < 0) CALL juDFT_error('not yet implemented', calledby='wannierlib_soc_tlmw')
      j = jwf(nwf)
      jm = jmwf(nwf)
      IF (j < 0) CALL juDFT_error('jwf', calledby='wannierlib_soc_tlmw')
      IF (ABS(jm) - j > 1e-10) CALL juDFT_error('jmwf', calledby='wannierlib_soc_tlmw')
      IF (ABS(l + 0.5 - j) > 1e-10 .AND. ABS(l - 0.5 - j) > 1e-10) CALL juDFT_error('regula trianguli violata', calledby='wannierlib_soc_tlmw')
      tlmwf(0:3, -3:3, nwf) = CMPLX(0.0, 0.0)
      DO m = -l, l
        tlmwf(l, m, nwf) = clebsch(REAL(l), 0.5, REAL(m), 1.5 - jspin, j, jm)
      END DO
    END DO
    CALL timestop('wannierlib_soc_tlmw')
  END SUBROUTINE wannierlib_soc_tlmw

END MODULE m_wannierlib_rad_twd
