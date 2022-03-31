!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_spinor
   IMPLICIT NONE
CONTAINS

   !The spinors are calculated both in hssphn_sph & hssphn_nonsph, hence this is a
   !common subroutine
  SUBROUTINE hsmt_spinor(iSpinNum, n, nococonv, chi_mat)
      USE m_types
      USE m_constants

      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      INTEGER,          INTENT(IN)  :: iSpinNum, n
      COMPLEX,          INTENT(OUT) :: chi_mat(2,2)

      INTEGER           :: iSpinPr, iSpin
      COMPLEX           :: umat(2,2)

      !--->       set up the spinors of this atom within global
      !--->       spin-coordinateframe
      umat = nococonv%umat(n)
      !--->       and determine the prefactors for the Hamitonian- and
      !--->       overlapp-matrix elements
      IF (iSpinNum<3) THEN
         iSpinPr = iSpinNum
         iSpin   = iSpinNum
      ELSE IF(iSpinNum==3) THEN
         iSpinPr = 2
         iSpin   = 1
      ELSE
         iSpinPr = 1
         iSpin   = 2
      ENDIF

      chi_mat(1, 1) = umat(1,iSpinPr)*CONJG(umat(1,iSpin))
      chi_mat(1, 2) = umat(1,iSpinPr)*CONJG(umat(2,iSpin))
      chi_mat(2, 1) = umat(2,iSpinPr)*CONJG(umat(1,iSpin))
      chi_mat(2, 2) = umat(2,iSpinPr)*CONJG(umat(2,iSpin))

   END SUBROUTINE hsmt_spinor

  SUBROUTINE hsmt_spinor_soc(n,ki,nococonv,lapw,chi_so,angso,kj_start,kj_end)
    USE m_types
    use m_constants

    TYPE(t_nococonv),INTENT(IN)      :: nococonv
    TYPE(t_lapw),INTENT(IN)      :: lapw
    INTEGER,INTENT(IN)           :: n,ki
    COMPLEX,INTENT(out)          :: chi_so(:,:,:,:)
    COMPLEX,INTENT(out),OPTIONAL :: angso(:,:,:)
    INTEGER,INTENT(in), OPTIONAL :: kj_start,kj_end

    REAL     :: cross_k(3)
    INTEGER  :: j1,j2,kj
    COMPLEX  :: isigma(2,2,3)
    COMPLEX  :: chi(2,2)
    COMPLEX  :: isigma_x(2,2),isigma_y(2,2),isigma_z(2,2),d(2,2)

    !     isigma= i * sigma, where sigma is Pauli matrix
    isigma=CMPLX(0.0,0.0)
    isigma(1,2,1)=CMPLX(0.0,1.0)
    isigma(2,1,1)=CMPLX(0.0,1.0)
    isigma(1,2,2)=CMPLX(-1.0,0.0)
    isigma(2,1,2)=CMPLX(1.0,0.0)
    isigma(1,1,3)=CMPLX(0.0,-1.0)
    isigma(2,2,3)=CMPLX(0.0,1.0)

    !--->       set up the spinors of this atom within global
    !--->       spin-coordinateframe
    chi=conjg(nococonv%umat(n))

    isigma_x=MATMUL(conjg(transpose(chi)), MATMUL(isigma(:,:,1),((chi))))
    isigma_y=MATMUL(conjg(transpose(chi)), MATMUL(isigma(:,:,2),((chi))))
    isigma_z=MATMUL(conjg(transpose(chi)), MATMUL(isigma(:,:,3),((chi))))
    DO j1=1,2
       DO j2=1,2
          chi_so(1,1,j1,j2)=chi(1,j1)*CONJG(chi(1,j2))
          chi_so(2,1,j1,j2)=chi(2,j1)*CONJG(chi(1,j2))
          chi_so(2,2,j1,j2)=chi(2,j1)*CONJG(chi(2,j2))
          chi_so(1,2,j1,j2)=chi(1,j1)*CONJG(chi(2,j2))
       ENDDO
    ENDDO
    IF (.not.present(angso)) RETURN !only chis are needed
    !In the first variation SOC case the off-diagonal spinors are needed
    IF (present(angso)) THEN
      IF ((.not.present(kj_start)).or.((.not.present(kj_end)))) RETURN
    ENDIF
    DO kj = kj_start,kj_end
       cross_k(1)=lapw%gk(2,ki,1)*lapw%gk(3,kj,1)- lapw%gk(3,ki,1)*lapw%gk(2,kj,1)
       cross_k(2)=lapw%gk(3,ki,1)*lapw%gk(1,kj,1)- lapw%gk(1,ki,1)*lapw%gk(3,kj,1)
       cross_k(3)=lapw%gk(1,ki,1)*lapw%gk(2,kj,1)- lapw%gk(2,ki,1)*lapw%gk(1,kj,1)
       DO j1=1,2
          DO j2=1,2
             angso(kj-kj_start+1,j1,j2)= (isigma_x(j1,j2)*cross_k(1)+&
                     isigma_y(j1,j2)*cross_k(2)+ isigma_z(j1,j2)*cross_k(3))
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE hsmt_spinor_soc


END MODULE m_hsmt_spinor
