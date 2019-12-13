!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_mtNocoPot_offdiag
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_mtNocoPot_offdiag(n,mpi,sym,atoms,noco,cell,lapw,td,fj,gj,hmat_tmp,hmat)
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
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_tlmplm),INTENT(IN)     :: td
#if defined CPP_GPU
    REAL,MANAGED,INTENT(IN)    :: fj(:,:,:,:),gj(:,:,:,:)
#else
    REAL,INTENT(IN)            :: fj(:,0:,:,:),gj(:,0:,:,:)
#endif
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)          :: n
    COMPLEX                       :: chi_one,chi(2,2)
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:),hmat_tmp

    chi_one=1.0
    CALL hmat_tmp%clear()
    !The spin1,2 matrix is calculated(real part of potential)
    CALL hsmt_nonsph(n,mpi,sym,atoms,3,1,1,chi_one,noco,cell,lapw,td,&
         fj(:,0:,1,:),gj(:,0:,1,:),hmat_tmp)

    CALL hsmt_spinor(3,n,noco,chi) !spinor for off-diagonal part
    CALL hsmt_distspins(chi,hmat_tmp,hmat)

    CALL hmat_tmp%TRANSPOSE()
    hmat_tmp%data_c=CONJG(hmat_tmp%data_c)
    CALL hsmt_spinor(4,n,noco,chi) !spinor for off-diagonal part
    CALL hsmt_distspins(chi,hmat_tmp,hmat)


    CALL hmat_tmp%clear()
    !The spin1,2 matrix is calculated(imag part of potential)
    chi_one=CMPLX(0.,1.)
    CALL hsmt_nonsph(n,mpi,sym,atoms,4,1,1,chi_one,noco,cell,lapw,td,&
         fj(:,0:,1,:),gj(:,0:,1,:),hmat_tmp)

    CALL hsmt_spinor(3,n,noco,chi) 
    CALL hsmt_distspins(chi,hmat_tmp,hmat)

    CALL hmat_tmp%TRANSPOSE()
    hmat_tmp%data_c=CONJG(hmat_tmp%data_c)
    CALL hsmt_spinor(4,n,noco,chi) 
    CALL hsmt_distspins(chi,hmat_tmp,hmat)
  END SUBROUTINE hsmt_mtNocoPot_offdiag
END MODULE m_hsmt_mtNocoPot_offdiag
