!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_elph_linewidth

#ifdef CPP_MPI
    USE mpi
#endif
    USE m_juDFT
    USE m_types 
    USE m_constants


    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_ph_linewidth(fi,fmpi,wtkpt,eig_k,eig_kq,gmat,eigenVals,ef,ph_linewidth)
        ! Phonon linewidth for one q-point via energy-grid binning.
        ! i_integration = 1 : single-delta (Allen transport) form with Fermi factors
        !               = 2 : double-delta approximation
        ! The eigenvalues, weights and matrix elements are passed in, so this works for
        ! both the coarse and the (Wannier-)interpolated meshes. The caller owns any file
        ! output. For a serial caller pass an fmpi whose mpi_comm is MPI_COMM_SELF, so the
        ! MPI_ALLREDUCE below is a no-op.

        USE m_dosbin
        USE m_smooth
        USE m_dfpt_fermie, ONLY : sfermi

        TYPE(t_fleurinput), INTENT(IN)  :: fi
        TYPE(t_mpi),        INTENT(IN)  :: fmpi
        REAL,               INTENT(IN)  :: wtkpt(:)          ! weight per k-point
        REAL,               INTENT(IN)  :: eig_k(:,:,:)      ! (nu,  nk, jsp) bands at k
        REAL,               INTENT(IN)  :: eig_kq(:,:,:)     ! (nu', nk, jsp) bands at k+q
        COMPLEX,            INTENT(IN)  :: gmat(:,:,:,:,:)   ! (nu',nu,k,jsp,mode); squared internally
        REAL,               INTENT(IN)  :: eigenVals(:)      ! omega^2 per mode
        REAL,               INTENT(IN)  :: ef
        REAL, ALLOCATABLE,  INTENT(OUT) :: ph_linewidth(:)   ! (nmode)

        INTEGER :: iNupr, nu, ispin, gridPoint, nk, nk_i, nZero, iMode, nmode, jspins, ndos
        REAL    :: emin, emax, x, xq
        REAL, ALLOCATABLE :: eGrid(:), linewidth(:,:)
        REAL, ALLOCATABLE :: gmat2(:,:,:,:,:)   ! |g|^2
#ifdef CPP_MPI
        INTEGER :: ierr
#endif

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
                    DO nk_i = 1, size(fmpi%k_list)
                        nk = fmpi%k_list(nk_i)
                        DO nu = 1, SIZE(gmat2,2)
                            x = (eig_k(nu,nk,ispin)-ef)/fi%input%tkb
                            DO iNupr = 1, SIZE(gmat2,1)
                                xq = (eig_kq(iNupr,nk,ispin)-ef)/fi%input%tkb
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
                    CALL dos_bin_transport(fmpi, jspins, wtkpt, eGrid, eig_k, eig_kq, &
                                           gmat2(:,:,:,:,iMode), linewidth, -SQRT(eigenVals(iMode)))
#ifdef CPP_MPI
                    CALL mpi_allreduce(MPI_IN_PLACE, linewidth, size(linewidth), MPI_DOUBLE_PRECISION, MPI_SUM, fmpi%mpi_comm, ierr)
#endif
                    DO ispin = 1, jspins
                        CALL smooth(eGrid, linewidth(:,ispin), fi%dfpt%smearingGauss, size(eGrid))
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
                    CALL dos_bin_double(fmpi, jspins, wtkpt, eGrid, eig_k, eig_kq, &
                                        gmat2(:,:,:,:,iMode), fi%dfpt%smearingGauss, linewidth, ef)
#ifdef CPP_MPI
                    CALL mpi_allreduce(MPI_IN_PLACE, linewidth, size(linewidth), MPI_DOUBLE_PRECISION, MPI_SUM, fmpi%mpi_comm, ierr)
#endif
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