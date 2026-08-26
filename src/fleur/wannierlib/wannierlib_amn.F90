!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  A_mn(k) = <psi_mk | g_n>, the overlap of each Bloch state with each trial orbital, which
!>  is what tells Wannier90 where to start.
!>
!>  Assembles what the other two halves produce: the radial functions from
!>  m_wannierlib_rad_twd, the angular coefficients from m_wannierlib_tlmw or its spin-orbit
!>  counterpart, and the augmentation coefficients of the state, which arrive already
!>  computed. Nothing here reads an eigenvector: the state reaches it as abc.
!>
!>  A trial orbital may be rotated, and the rotation is applied to its angular coefficients
!>  with Wigner D matrices rather than to the harmonics themselves.
!>
!>  Fills one k-point's column of amn, which the caller sums over its own k-slice.
MODULE m_wannierlib_amn
  USE m_juDFT
  USE m_constants
  USE m_intgr, ONLY: intgr3
  USE m_dwigner
  USE m_eulerrot
  USE m_clebsch
  USE m_types_atoms
  USE m_types_kpts
  USE m_types_abc
  USE m_types_usdus
  USE m_types_radfun
  USE m_types_wannierlib
  USE m_wannierlib_tlmw
  USE m_wannierlib_rad_twd
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_amn
CONTAINS

  SUBROUTINE wannierlib_amn(wannierlib, atoms, kpts, ikpt, usdus, radfun, abc, l_nocosoc, l_spinors, jspin, jspin_rad, amn)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: wannierlib
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_kpts), INTENT(IN) :: kpts
    INTEGER, INTENT(IN) :: ikpt
    TYPE(t_usdus), INTENT(IN) :: usdus
    TYPE(t_radfun), INTENT(IN) :: radfun(atoms%ntype)
    TYPE(t_abc), INTENT(IN) :: abc(atoms%ntype)
    LOGICAL, INTENT(IN) :: l_nocosoc
    LOGICAL, INTENT(IN) :: l_spinors
    INTEGER, INTENT(IN) :: jspin       ! physical spin (filters the projections)
    INTEGER, INTENT(IN) :: jspin_rad   ! radial index (=1 when jspins=1)
    COMPLEX, INTENT(INOUT) :: amn(:,:)
    
    INTEGER :: nwf, nat_local, ntyp, ne, l, m, lm, j, mp, ir
    REAL :: rads(wannierlib%num_wann, 0:3, atoms%jmtd, 2)
    REAL :: vlpr(atoms%jmtd)
    REAL :: proj_int(2 + atoms%nlod)
    COMPLEX :: tlmwf(0:3, -3:3, wannierlib%num_wann), tlmwft(0:3, -3:3, wannierlib%num_wann)
    REAL :: amx(3, 3, wannierlib%num_wann), imx(3, 3)
    COMPLEX :: wign(-3:3, -3:3, 3, wannierlib%num_wann)
    COMPLEX :: factor, coeff_sum
    REAL :: arg
    LOGICAL :: has_soc_proj

    IF (wannierlib%num_wann <= 0) THEN
      CALL juDFT_error('wannierlib_amn: no projections configured', calledby='wannierlib_amn')
    END IF
    
    CALL timestart('wannierlib_amn')

    CALL wannierlib_rad_twd(wannierlib, atoms, wannierlib%num_wann, ikpt, usdus, radfun, jspin_rad, rads)

    tlmwf = CMPLX(0.0, 0.0)
    tlmwft = CMPLX(0.0, 0.0)

    has_soc_proj = ALL(wannierlib%proj_j(1:wannierlib%num_wann) > 0.0)
    IF (l_nocosoc .AND. has_soc_proj) THEN
      CALL wannierlib_soc_tlmw(wannierlib%num_wann, wannierlib%proj_l, wannierlib%proj_j, wannierlib%proj_mj, jspin, tlmwf)
    ELSE
      CALL wannierlib_tlmw(wannierlib, wannierlib%num_wann, l_nocosoc, l_spinors, jspin, tlmwf)
    END IF

    CALL eulerrot(wannierlib%num_wann, wannierlib%proj_alpha, wannierlib%proj_beta, wannierlib%proj_gamma, amx)
    imx(:, :) = 0.0
    imx(1, 1) = 1.0
    imx(2, 2) = 1.0
    imx(3, 3) = 1.0
    CALL d_wigner(wannierlib%num_wann, amx, imx, 3, wign)

    DO nwf = 1, wannierlib%num_wann
      tlmwft(0, :, nwf) = tlmwf(0, :, nwf)
      DO l = 1, 3
        DO m = -l, l
          DO mp = -l, l
            tlmwft(l, m, nwf) = tlmwft(l, m, nwf) + wign(mp, m, l, nwf) * tlmwf(l, mp, nwf)
          END DO
        END DO
      END DO
    END DO
    
    DO nwf = 1, wannierlib%num_wann
      arg = -kpts%bkf(1, ikpt) * wannierlib%proj_shift(1, nwf)
      arg = arg - kpts%bkf(2, ikpt) * wannierlib%proj_shift(2, nwf)
      arg = arg - kpts%bkf(3, ikpt) * wannierlib%proj_shift(3, nwf)
      arg = tpi_const * arg
      factor = CMPLX(COS(arg), SIN(arg)) * wannierlib%proj_weight(nwf)
    

      ntyp = wannierlib%proj_ntype(nwf)
      nat_local = wannierlib%proj_atom(nwf)
      DO ne = 1, size(amn,1)
        DO l = 0, MIN(atoms%lmax(ntyp), 3)
          proj_int(:) = 0.0
          DO j = 1, abc(ntyp)%n_r(l)
            DO ir = 1, atoms%jri(ntyp)
              vlpr(ir) = radfun(ntyp)%r(ir, 1, j, l, jspin_rad) * rads(nwf, l, ir, 1) + &
                         radfun(ntyp)%r(ir, 2, j, l, jspin_rad) * rads(nwf, l, ir, 2)
              IF (wannierlib%proj_rwf(nwf) > 0) THEN
                vlpr(ir) = vlpr(ir) * atoms%rmsh(ir, ntyp)
              END IF
            END DO
            CALL intgr3(vlpr, atoms%rmsh(1, ntyp), atoms%dx(ntyp), atoms%jri(ntyp), proj_int(j))
          END DO

          DO m = -l, l
            lm = l * (l + 1) + m
            coeff_sum = CMPLX(0.0, 0.0)
            DO j = 1, abc(ntyp)%n_r(l)
              coeff_sum = coeff_sum + abc(ntyp)%cof(ne, lm, j, nat_local) * proj_int(j)
            END DO
            amn(ne, nwf) = amn(ne, nwf) + tlmwft(l, m, nwf) * CONJG(coeff_sum * (ImagUnit) ** l) * factor
          END DO
        END DO
      END DO
    END DO

    CALL timestop('wannierlib_amn')
  END SUBROUTINE wannierlib_amn

 

END MODULE m_wannierlib_amn
