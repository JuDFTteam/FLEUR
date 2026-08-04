!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Provider of the orbital angular momentum operator L in the Bloch basis,
!>      L0_alpha,mn(k) = < psi_mk | L_alpha | psi_nk > ,   alpha = x,y,z.
!>
!>  Class-A (on-site) operator: evaluated inside the muffin-tin spheres from the
!>  angular-momentum matrix of the complex spherical harmonics (mirrors FLEUR's
!>  wann_anglmom): in the Ylm basis L_z is diagonal (= m), and L_+/L_- connect
!>  m -> m+1 / m -> m-1 with sqrt((l-+m)(l+-m+1)). L is orbital and spin-DIAGONAL,
!>  so we simply sum both local spin components abc(:,1)+abc(:,2) -- no ccchi
!>  rotation is needed (the spin trace is frame-invariant). The radial part is the
!>  ordinary overlap radfun%integral (L does not touch it). Interstitial and vacuum
!>  contributions to L are neglected (orbital moment is a muffin-tin quantity).
MODULE m_melem_orbmom
  USE m_juDFT
  USE m_constants, ONLY : ImagUnit
  USE m_types_atoms
  USE m_types_abc
  USE m_types_radfun
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_orbmom_bloch_collinear
CONTAINS

  !> Site-summed L0(:,:,1:3) at one k for a single spin channel: with jspins=2 and no
  !> SOC or noco the two channels are separate scalar problems, so the coefficients
  !> arrive for one channel only (abc(:), radial index jspin_rad) and there is nothing
  !> to sum over spin. Muffin-tin only, as the orbital moment is.
  SUBROUTINE melem_orbmom_bloch_collinear(atoms, abc, radfun, jspin_rad, l0)
    TYPE(t_atoms),  INTENT(IN)  :: atoms
    TYPE(t_abc),    INTENT(IN)  :: abc(:)         ! (ntype) single-channel local-frame coeffs
    TYPE(t_radfun), INTENT(IN)  :: radfun(:)      ! (ntype) : %integral(.,.,l,jspins,jspins)
    INTEGER,        INTENT(IN)  :: jspin_rad      ! radial spin slot of this channel (<= jspins)
    COMPLEX,        INTENT(OUT) :: l0(:, :, :)    ! (nb,nb,3): 1=Lx 2=Ly 3=Lz (site-summed)

    INTEGER :: nb, i, j, ntyp, iat, l, ll1, mm, lm, n_r, n_r2, s
    REAL    :: lplus, lminus, w
    COMPLEX :: cz, cp, cm    ! L_z, L_+, L_- accumulators for (i,j), summed over atoms

    nb = SIZE(l0, 1); s = jspin_rad
    l0 = CMPLX(0.0, 0.0)
    DO j = 1, nb                       ! ket band
      DO i = 1, nb                     ! bra band
        cz = CMPLX(0.0, 0.0); cp = CMPLX(0.0, 0.0); cm = CMPLX(0.0, 0.0)
        DO ntyp = 1, atoms%ntype
          DO iat = 1, atoms%neq(ntyp)
            DO l = 0, atoms%lmax(ntyp)
              ll1 = l*(l + 1)
              DO mm = -l, l
                lm = ll1 + mm
                lplus  = SQRT(REAL((l - mm)*(l + mm + 1)))   ! <m+1|L+|m>
                lminus = SQRT(REAL((l + mm)*(l - mm + 1)))   ! <m-1|L-|m>
                DO n_r = 1, abc(ntyp)%n_r(l)
                  DO n_r2 = 1, abc(ntyp)%n_r(l)
                    w = radfun(ntyp)%integral(n_r, n_r2, l, s, s)
                    cz = cz + abc(ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(ntyp)%cof(j,lm,n_r2,iat))*REAL(mm)*w
                    IF (mm < l) cp = cp + abc(ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(ntyp)%cof(j,lm+1,n_r2,iat))*lplus*w
                    IF (mm > -l) cm = cm + abc(ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(ntyp)%cof(j,lm-1,n_r2,iat))*lminus*w
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        l0(i, j, 1) = 0.5 * (cp + cm)              ! Lx = (L+ + L-)/2
        l0(i, j, 2) = -0.5 * ImagUnit * (cp - cm)  ! Ly = (L+ - L-)/(2i)
        l0(i, j, 3) = cz                           ! Lz
      END DO
    END DO
  END SUBROUTINE melem_orbmom_bloch_collinear

END MODULE m_melem_orbmom
