!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Provider of the spin-orbit-coupling operator in the Bloch basis,
!>      H0_SOC,mn(k) = < psi_mk | xi(r) L.S | psi_nk > .
!>
!>  Class-A (on-site) operator, but unlike spin/L it needs the RELATIVISTIC radial
!>  SOC integrals xi(r)=(1/2m^2c^2)(1/r)dV/dr. Rather than re-derive them we reuse
!>  FLEUR's dedicated machinery (adapts wann_socmat): spnorb builds the radial SOC
!>  integrals + the L.S angular/spin matrix (rsoc), and hsoham contracts them with
!>  the wavefunction coefficients into the SOC Hamiltonian block hsomtx(nb,nb,2,2).
!>  The full spinor matrix element is the sum over the 2x2 spin blocks.
!>
!>  hsoham does MPI_BARRIER/MPI_REDUCE on SUB_COMM; since this runs on the master
!>  rank only we pass MPI_COMM_SELF (n_size=1, all atoms local) so the collectives
!>  are trivial. Local orbitals are neglected in this first version (chelp = 0).
MODULE m_wannierlib_socmat_melem
  USE m_juDFT
  USE m_spnorb
  USE m_hsoham
  USE m_types
  USE m_types_abc
#ifdef CPP_MPI
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_socmat_bloch
  ! spnorb is k-independent: compute rsoc once and cache it (also avoids 512x printing)
  TYPE(t_rsoc), SAVE :: rsoc_cache
  LOGICAL,      SAVE :: rsoc_ready = .FALSE.
CONTAINS

  !> Build H0_SOC(:,:) (single component) in the Bloch basis at one k.
  SUBROUTINE wannierlib_socmat_bloch(atoms, noco, nococonv, input, fmpi, enpara, vtot, usdus, abc, nb, soc0, soc4)
    TYPE(t_atoms),    INTENT(IN)    :: atoms
    TYPE(t_noco),     INTENT(IN)    :: noco
    TYPE(t_nococonv), INTENT(IN)    :: nococonv
    TYPE(t_input),    INTENT(IN)    :: input
    TYPE(t_mpi),      INTENT(IN)    :: fmpi
    TYPE(t_enpara),   INTENT(IN)    :: enpara
    TYPE(t_potden),   INTENT(IN)    :: vtot
    TYPE(t_usdus),    INTENT(INOUT) :: usdus
    TYPE(t_abc),      INTENT(IN)    :: abc(:, :)     ! (ntype, 2 spin) local-frame coeffs
    INTEGER,          INTENT(IN)    :: nb
    COMPLEX,          INTENT(OUT)   :: soc0(:, :, :) ! (nb,nb,1) SOC operator (single component)
    COMPLEX,          INTENT(OUT)   :: soc4(:, :, :) ! (nb,nb,4) 2x2 SOC spin blocks, c=(ii-1)*2+jj (rssocmat.1)

    COMPLEX, ALLOCATABLE :: ahelp(:, :, :, :), bhelp(:, :, :, :), chelp(:, :, :, :, :)
    COMPLEX, ALLOCATABLE :: hsomtx(:, :, :, :)
    INTEGER :: lmd, natd, nsz(2), ntyp, iat, na, l, ll1, m, lm, ie, isp, comm_self

    lmd  = atoms%lmaxd*(atoms%lmaxd + 2)
    natd = atoms%nat
    nsz  = nb

    ! ---- pack a/b coefficients as ahelp(lm,na,ie,ispin); LO neglected (chelp=0) ----
    ALLOCATE(ahelp(lmd, natd, nb, input%jspins), source=CMPLX(0.0, 0.0))
    ALLOCATE(bhelp(lmd, natd, nb, input%jspins), source=CMPLX(0.0, 0.0))
    ALLOCATE(chelp(-atoms%llod:atoms%llod, nb, atoms%nlod, natd, input%jspins), source=CMPLX(0.0, 0.0))
    DO isp = 1, input%jspins
      na = 0
      DO ntyp = 1, atoms%ntype
        DO iat = 1, atoms%neq(ntyp)
          na = na + 1
          DO l = 1, atoms%lmax(ntyp)         ! SOC: l >= 1
            ll1 = l*(l + 1)
            DO m = -l, l
              lm = ll1 + m
              DO ie = 1, nb
                ahelp(lm, na, ie, isp) = abc(ntyp, isp)%cof(ie, lm, 1, iat)   ! u
                bhelp(lm, na, ie, isp) = abc(ntyp, isp)%cof(ie, lm, 2, iat)   ! udot
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ! ---- relativistic radial SOC integrals + L.S angular/spin matrix (computed ONCE) ----
    IF (.NOT. rsoc_ready) THEN
      ! spnorb (via sorad) OVERWRITES usdus. Use a PRIVATE copy so the shared usdus that
      ! feeds the wannierization (amn/mmn) is left intact -> reproducible spread with SOC.
      BLOCK
        TYPE(t_usdus) :: usdus_soc
        usdus_soc = usdus
        CALL spnorb(atoms, noco, nococonv, input, fmpi, enpara, vtot%mt, usdus_soc, rsoc_cache, .TRUE.)
      END BLOCK
      rsoc_cache%soangl = CONJG(rsoc_cache%soangl)
      rsoc_ready = .TRUE.
    END IF

    ! ---- SOC Hamiltonian block hsomtx(nb,nb,2,2); master-rank-only via COMM_SELF ----
    comm_self = fmpi%mpi_comm
#ifdef CPP_MPI
    comm_self = MPI_COMM_SELF
#endif
    ALLOCATE(hsomtx(nb, nb, 2, 2), source=CMPLX(0.0, 0.0))
    CALL hsoham(atoms, noco, input, nsz, nb, chelp, rsoc_cache, ahelp, bhelp, &
                1, natd, 0, 1, comm_self, hsomtx)

    ! ---- full spinor matrix element = sum over the 2x2 spin blocks ----
    soc0(:, :, 1) = hsomtx(:, :, 1, 1) + hsomtx(:, :, 1, 2) + hsomtx(:, :, 2, 1) + hsomtx(:, :, 2, 2)
    ! the 4 spin blocks (uncollapsed) for the real-space SOC export rssocmat.1: c=(ii-1)*2+jj
    soc4(:, :, 1) = hsomtx(:, :, 1, 1)   ! jj=1, ii=1
    soc4(:, :, 2) = hsomtx(:, :, 2, 1)   ! jj=2, ii=1
    soc4(:, :, 3) = hsomtx(:, :, 1, 2)   ! jj=1, ii=2
    soc4(:, :, 4) = hsomtx(:, :, 2, 2)   ! jj=2, ii=2

    DEALLOCATE(ahelp, bhelp, chelp, hsomtx)
  END SUBROUTINE wannierlib_socmat_bloch

END MODULE m_wannierlib_socmat_melem
