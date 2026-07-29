!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Spin-orbit coupling matrix element in the Bloch basis,
!>
!>      SOC_mn(k) = < psi_mk | H_SOC | psi_nk > ,
!>
!>  as the single collapsed component (soc0) and as the four uncollapsed 2x2 spinor
!>  blocks (soc4, written to rssocmat.1).
!>
!>  !! CURRENTLY UNAVAILABLE -- awaiting a refactor onto the new SOC infrastructure. !!
!>
!>  This provider used to build the operator from m_spnorb (relativistic radial SOC
!>  integrals + the L.S angular matrix soangl) and m_hsoham (the Bloch-basis SOC
!>  Hamiltonian). Both routines were removed with the SOC refactor that replaced
!>  alineso/eigenso by m_secvar_soc, so the dependency is dropped here rather than
!>  carried on removed code.
!>
!>  The replacement path is m_types_matelements_soc: t_matelements_soc%calc_matrix_elements
!>  builds exactly these 2x2 spinor blocks from the abc coefficients and a t_rsoc
!>  (t_rsoc%init / %rad_matrix / %angles), and it even includes the local-orbital
!>  channels that the old m_hsoham call neglected. Porting to it requires
!>    * the abc coefficients in (2,ntype) order (this layer holds them as (ntype,2)),
!>    * a t_rsoc built once per potential rather than the previous cached spnorb call,
!>    * settling whether the CONJG(soangl) fix-up that this routine used to apply is
!>      still needed -- m_hsoham and calc_matrix_elements share the
!>      soangl * conjg(cof) convention, so that has to be re-derived and validated
!>      against the wannierlib SOC tests, not guessed.
!>  Until that is done the operator reports itself as unavailable instead of
!>  returning silently wrong numbers.
MODULE m_melem_socmat
  USE m_juDFT
  USE m_types
  USE m_types_abc
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_socmat_bloch
  LOGICAL, SAVE :: warned = .FALSE.
CONTAINS

  !> Build H0_SOC(:,:) in the Bloch basis at one k. See the module header: the
  !> implementation is pending the port onto m_types_matelements_soc, so this returns
  !> zeros and warns once. The argument list is the one the future implementation needs,
  !> so callers do not have to change when it lands.
  SUBROUTINE melem_socmat_bloch(atoms, noco, nococonv, input, fmpi, enpara, vtot, usdus, abc, nb, soc0, soc4)
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
    COMPLEX,          INTENT(OUT)   :: soc4(:, :, :) ! (nb,nb,4) 2x2 SOC spin blocks, c=(ii-1)*2+jj

    soc0 = CMPLX(0.0, 0.0)
    soc4 = CMPLX(0.0, 0.0)

    IF (.NOT. warned) THEN
      warned = .TRUE.
      CALL juDFT_warn("wannierlib: the SOC matrix element is not available in this build", &
                      hint="m_melem_socmat still has to be ported onto m_types_matelements_soc "// &
                           "after the SOC refactor removed m_spnorb/m_hsoham. Any requested "// &
                           "'soc'/'spin_orbit' operator output (bands_wann_soc.dat, rssocmat.1) is zero.", &
                      calledby="melem_socmat_bloch")
    END IF
  END SUBROUTINE melem_socmat_bloch

END MODULE m_melem_socmat
