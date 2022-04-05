!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_mtNocoPot_offdiag
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_mtNocoPot_offdiag(n,input,fmpi,sym,atoms,noco,nococonv,cell,lapw,ud,td,fjgj,igSpinPr,igSpin,hmat_tmp,hmat)
    !Calculate the contribution from the local-spin-offdiagonal potential
    !The following idea is used:
    !Calculate the matrix by using non-spherical algorithm. This is done only once, since
    !this sets up both the local spin-up-down and the spin-down-up part (it calculates the
    !full matrix). So both can be updated from this matrix. But since the off-diagonal
    !local potential is real we have to call the routine twice and use the chi_one factor
    !to get the imaginary contribution
    USE m_types
    USE m_hsmt_nonsph
    USE m_hsmt_distspins
    USE m_hsmt_spinor
    USE m_hsmt_lo
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_nococonv),INTENT(IN)       :: nococonv
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_usdus),INTENT(IN)      :: ud
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_fjgj),INTENT(IN)       :: fjgj
    INTEGER,INTENT(IN)            :: igSpinPr,igSpin

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)          :: n
    COMPLEX                       :: chi_one,chi(2,2)
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:),hmat_tmp

    chi_one=1.0
    !The spin2,1 matrix is calculated(real part of potential)
    CALL hsmt_nonsph(n,fmpi,sym,atoms,2,1,igSpinPr,igSpin,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,.FALSE.,1)
    CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,ud,td,fjgj,n,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.TRUE.,.FALSE.,.FALSE.)
    !call hmat_tmp%u2l()
    CALL hsmt_spinor(3,n,nococonv,chi) !spinor for off-diagonal part
    CALL hsmt_distspins(chi,hmat_tmp,hmat)

    !hmat_tmp%data_c=conjg(hmat_tmp%data_c)
    !CALL hmat_tmp%TRANSPOSE()
    !CALL hsmt_spinor(4,n,nococonv,chi) !spinor for off-diagonal part
    !CALL hsmt_distspins(chi,hmat_tmp,hmat)


    !The spin1,2 matrix is calculated(imag part of potential)
    !chi_one=CMPLX(0.,1.)
    CALL hsmt_nonsph(n,fmpi,sym,atoms,1,2,igSpinPr,igSpin,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,.FALSE.,1)
    CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,ud,td,fjgj,n,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.TRUE.,.FALSE.,.FALSE.)
    !call hmat_tmp%u2l()

    CALL hsmt_spinor(4,n,nococonv,chi)
    CALL hsmt_distspins(chi,hmat_tmp,hmat)
    !hmat_tmp%data_c=conjg(hmat_tmp%data_c)
    !CALL hmat_tmp%TRANSPOSE()
    !CALL hsmt_spinor(4,n,nococonv,chi)
    !CALL hsmt_distspins(chi,hmat_tmp,hmat)
  END SUBROUTINE hsmt_mtNocoPot_offdiag
END MODULE m_hsmt_mtNocoPot_offdiag
