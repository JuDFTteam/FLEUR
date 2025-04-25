!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_spinor
   !!Module for the calculation of transformation matrices between 
   !!Spin-dependent matrices in the local frame and the global frame
   IMPLICIT NONE
CONTAINS


  SUBROUTINE hsmt_spinor(iSpinNum, n, nococonv, chi_mat)
   !! Obtain the part of the transformation that maps the local spin
   !! iSpinNum to the global spin frame. 
      USE m_types
      USE m_constants

      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      INTEGER,          INTENT(IN)  :: iSpinNum     !! local spin
      INTEGER,          INTENT(IN)  :: n            !! atom number
      COMPLEX,          INTENT(OUT) :: chi_mat(2,2) !! Matrix describing the mapping

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
      ELSE IF(iSpinNum==4) THEN
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

  SUBROUTINE hsmt_spinor_soc(n,nococonv,chi_so,isigma_xyz)
    !$acc routine seq
    !!Generalization of hsmt_spinor to SOC case. 
    USE m_types
    use m_constants

    TYPE(t_nococonv),INTENT(IN)  :: nococonv
    INTEGER,INTENT(IN)           :: n   !!index of atom
    COMPLEX,INTENT(out)          :: chi_so(:,:,:,:) !! Transformation from local to global spin
    COMPLEX,INTENT(OUT),optional :: isigma_xyz(:,:,:)
    
    INTEGER  :: j1,j2
    COMPLEX  :: isigma(2,2,3)
    COMPLEX  :: chi(2,2)


    !--->       set up the spinors of this atom within global
    !--->       spin-coordinateframe
    chi=chi_explicit_nopass(nococonv%alph(n),nococonv%beta(n))
    
    if (present(isigma_xyz)) THEN
      !     isigma= i * sigma, where sigma is Pauli matrix
      isigma(:,:,:) = CMPLX(0.0,0.0)

      isigma(1,2,1)=CMPLX(0.0,1.0)  !     (0  1)   ( 0  i)
      isigma(2,1,1)=CMPLX(0.0,1.0)  ! i * (1  0) = ( i  0)
      isigma(1,2,2)=CMPLX(1.0,0.0)  !     (0 -i)   ( 0  1) 
      isigma(2,1,2)=CMPLX(-1.0,0.0) ! i * (i  0) = (-1  0) 
      isigma(1,1,3)=CMPLX(0.0,1.0)  !     (1  0)   ( i  0)
      isigma(2,2,3)=CMPLX(0.0,-1.0) ! i * (0 -1) = ( 0 -i)

      isigma_xyz(:,:,1)=MATMUL(conjg(transpose(chi)), MATMUL(isigma(:,:,1),chi))
      isigma_xyz(:,:,2)=MATMUL(conjg(transpose(chi)), MATMUL(isigma(:,:,2),chi))
      isigma_xyz(:,:,3)=MATMUL(conjg(transpose(chi)), MATMUL(isigma(:,:,3),chi))
    ENDIF  
   
    !chi=conjg(nococonv%umat(n))
    DO j1=1,2
       DO j2=1,2
          chi_so(1,1,j1,j2)=chi(1,j1)*CONJG(chi(1,j2))
          chi_so(2,1,j1,j2)=chi(2,j1)*CONJG(chi(1,j2))
          chi_so(2,2,j1,j2)=chi(2,j1)*CONJG(chi(2,j2))
          chi_so(1,2,j1,j2)=chi(1,j1)*CONJG(chi(2,j2))
       ENDDO
    ENDDO
    

  END SUBROUTINE hsmt_spinor_soc


END MODULE m_hsmt_spinor
