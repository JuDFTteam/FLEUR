!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_elph_linewidth

    USE m_juDFT
    USE m_types
    USE m_constants


    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_ph_linewidth(fi,wtkpt,eig_k,eig_kq,gmat,eigenVals,ef,ph_linewidth)
        ! Phonon linewidth for one q-point via energy-grid binning.
        ! i_integration = 1 : single-delta (Allen transport) form with Fermi factors
        !               = 2 : double-delta approximation
        ! All arrays are indexed with the LOCAL k index (1..size(gmat,3))

        USE m_dosbin
        USE m_smooth
        USE m_dfpt_fermie, ONLY : sfermi

        TYPE(t_fleurinput), INTENT(IN)  :: fi
        REAL,               INTENT(IN)  :: wtkpt(:)          ! weight per k-point
        REAL,               INTENT(IN)  :: eig_k(:,:,:)      ! (nu,  nk, jsp) bands at k
        REAL,               INTENT(IN)  :: eig_kq(:,:,:)     ! (nu', nk, jsp) bands at k+q
        COMPLEX,            INTENT(IN)  :: gmat(:,:,:,:,:)   ! (nu',nu,k,jsp,mode); squared internally
        REAL,               INTENT(IN)  :: eigenVals(:)      ! omega^2 per mode
        REAL,               INTENT(IN)  :: ef
        REAL, ALLOCATABLE,  INTENT(OUT) :: ph_linewidth(:)   ! (nmode)

        INTEGER :: iNupr, nu, ispin, gridPoint, nk_i, nZero, iMode, nmode, jspins, ndos
        REAL    :: emin, emax, x, xq
        REAL, ALLOCATABLE :: eGrid(:), linewidth(:,:)
        REAL, ALLOCATABLE :: gmat2(:,:,:,:,:)   ! |g|^2

        jspins = fi%input%jspins
        nmode  = SIZE(eigenVals)
        ndos   = fi%banddos%ndos_points

        ! |g|^2 (do not mutate the caller's gmat)
        ALLOCATE(gmat2(SIZE(gmat,1),SIZE(gmat,2),SIZE(gmat,3),SIZE(gmat,4),SIZE(gmat,5)))
        gmat2 = ABS(gmat)**2

        emin = -8 * fi%input%tkb
        emax =  8 * fi%input%tkb
        ALLOCATE(eGrid(ndos), linewidth(ndos,jspins), ph_linewidth(nmode))
        ph_linewidth = 0.0
        DO gridPoint = 1, ndos
            eGrid(gridPoint) = emin + (emax-emin)/(ndos-1.0)*(gridPoint-1.0)
        END DO
        nZero = MINLOC(ABS(eGrid), dim=1)   ! grid point closest to E_F

        SELECT CASE(fi%dfpt%i_integration)

            CASE(1)   ! single delta with Fermi occupation factors
                ! multiply |g|^2 by (f(eps_k) - f(eps_k+q)) so only occ -> unocc scattering contributes
                CALL timestart("Fermi factor multiplication")
                DO ispin = 1, jspins
                    DO nk_i = 1, SIZE(gmat2,3)
                        DO nu = 1, SIZE(gmat2,2)
                            x = (eig_k(nu,nk_i,ispin)-ef)/fi%input%tkb
                            DO iNupr = 1, SIZE(gmat2,1)
                                xq = (eig_kq(iNupr,nk_i,ispin)-ef)/fi%input%tkb
                                gmat2(iNupr,nu,nk_i,ispin,:) = gmat2(iNupr,nu,nk_i,ispin,:)*(sfermi(x) - sfermi(xq))
                            END DO ! iNupr
                        END DO ! nu
                    END DO ! nk
                END DO ! ispin
                CALL timestop("Fermi factor multiplication")

                CALL timestart("Single-delta binning")
                DO iMode = 1, nmode
                    IF (eigenVals(iMode) < 0.0) CYCLE   ! imaginary mode -> linewidth 0
                    linewidth = 0.0
                    ! delta(eps_k - eps_k+q - omega); omega = sqrt(eigenVals)
                    CALL dos_bin_transport(jspins, wtkpt, eGrid, eig_k, eig_kq, &
                                           gmat2(:,:,:,:,iMode), linewidth, -SQRT(eigenVals(iMode)))
                    DO ispin = 1, jspins
                        CALL smooth(eGrid, linewidth(:,ispin), fi%input%tkb, size(eGrid))
                        ph_linewidth(iMode) = ph_linewidth(iMode) + tpi_const*linewidth(nZero,ispin)
                    END DO
                END DO
                CALL timestop("Single-delta binning")

            CASE(2)   ! double delta
                CALL timestart("Double-delta binning")
                DO iMode = 1, nmode
                    IF (eigenVals(iMode) < 0.0) CYCLE   ! imaginary mode -> linewidth 0
                    linewidth = 0.0
                    ! delta(eps_k - E_F) * delta(eps_k+q - E_F); Gaussian smearing done inside
                    CALL dos_bin_double(jspins, wtkpt, eGrid, eig_k, eig_kq, &
                                        gmat2(:,:,:,:,iMode), fi%input%tkb, linewidth, ef)
                    DO ispin = 1, jspins
                        ! factor two for spin deg. is calculated in dos_bin
                        ph_linewidth(iMode) = ph_linewidth(iMode) + tpi_const*SQRT(eigenVals(iMode))*linewidth(nZero,ispin)
                    END DO
                END DO
                CALL timestop("Double-delta binning")

            CASE DEFAULT
                CALL juDFT_error("dfpt_ph_linewidth: only i_integration=1 (single delta) or 2 (double delta) supported", &
                                 calledby="dfpt_ph_linewidth")

        END SELECT

    END SUBROUTINE dfpt_ph_linewidth




END MODULE m_dfpt_elph_linewidth