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
  SUBROUTINE hsmt_spinor(isp,n,noco,chi_mat)
    USE m_types
    IMPLICIT NONE

    TYPE(t_noco),INTENT(IN)      :: noco
    INTEGER,INTENT(IN)           :: isp, n
    COMPLEX,INTENT(OUT)          :: chi_mat(2,2)
   
    INTEGER           :: isp1,isp2
    COMPLEX,PARAMETER :: ci=CMPLX(0.0,1.0)
    COMPLEX           :: chi(2,2)

    !--->       set up the spinors of this atom within global
    !--->       spin-coordinateframe
    chi(1,1) =  exp(-ci*noco%alph(n)/2)*cos(noco%beta(n)/2)
    chi(1,2) = -exp(-ci*noco%alph(n)/2)*sin(noco%beta(n)/2)
    chi(2,1) =  exp(ci*noco%alph(n)/2)*sin(noco%beta(n)/2)
    chi(2,2) =  exp(ci*noco%alph(n)/2)*cos(noco%beta(n)/2)
    !--->       and determine the prefactors for the Hamitonian- and
    !--->       overlapp-matrix elements
    chi_mat(1,1) = chi(1,isp)*CONJG(chi(1,isp))
    chi_mat(2,1) = chi(2,isp)*CONJG(chi(1,isp))
    chi_mat(2,2) = chi(2,isp)*CONJG(chi(2,isp))
    chi_mat(1,2) = chi(1,isp)*CONJG(chi(2,isp))

    
    
  END SUBROUTINE hsmt_spinor

  SUBROUTINE hsmt_spinor_soc(n,ki,noco,lapw,chi_so,angso)
    USE m_types
    IMPLICIT NONE
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_lapw),INTENT(IN)      :: lapw
    INTEGER,INTENT(IN)           :: n,ki
    COMPLEX,INTENT(out)          :: chi_so(:,:,:,:)
    COMPLEX,INTENT(out),OPTIONAL :: angso(:,:,:)

    REAL     :: cross_k(3)
    INTEGER  :: j1,j2,kj
    COMPLEX  :: isigma(2,2,3)
    COMPLEX  :: chi(2,2)
    COMPLEX  :: isigma_x(2,2),isigma_y(2,2),isigma_z(2,2)
    COMPLEX,PARAMETER :: ci=CMPLX(0.0,1.0)
    
    !     isigma= -i * sigma, where sigma is Pauli matrix
    isigma=CMPLX(0.0,0.0)
    isigma(1,2,1)=CMPLX(0.0,-1.0)
    isigma(2,1,1)=CMPLX(0.0,-1.0)
    isigma(1,2,2)=CMPLX(-1.0,0.0)
    isigma(2,1,2)=CMPLX(1.0,0.0)
    isigma(1,1,3)=CMPLX(0.0,-1.0)
    isigma(2,2,3)=CMPLX(0.0,1.0)
    
    !--->       set up the spinors of this atom within global
    !--->       spin-coordinateframe
    chi(1,1) =  exp(-ci*noco%alph(n)/2)*cos(noco%beta(n)/2)
    chi(1,2) = -exp(-ci*noco%alph(n)/2)*sin(noco%beta(n)/2)
    chi(2,1) =  exp(ci*noco%alph(n)/2)*sin(noco%beta(n)/2)
    chi(2,2) =  EXP(ci*noco%alph(n)/2)*COS(noco%beta(n)/2)

    isigma_x=MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(isigma(:,:,1),chi))
    isigma_y=MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(isigma(:,:,2),chi))
    isigma_z=MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(isigma(:,:,3),chi))
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
       DO kj = 1,ki
          cross_k(1)=lapw%gk(2,ki,1)*lapw%gk(3,kj,1)- lapw%gk(3,ki,1)*lapw%gk(2,kj,1)
          cross_k(2)=lapw%gk(3,ki,1)*lapw%gk(1,kj,1)- lapw%gk(1,ki,1)*lapw%gk(3,kj,1)
          cross_k(3)=lapw%gk(1,ki,1)*lapw%gk(2,kj,1)- lapw%gk(2,ki,1)*lapw%gk(1,kj,1)
          DO j1=1,2
             DO j2=1,2
                angso(kj,j1,j2)= isigma_x(j1,j2)*cross_k(1)+&
                     isigma_y(j1,j2)*cross_k(2)+ isigma_z(j1,j2)*cross_k(3)
             ENDDO
          ENDDO
       ENDDO
     END SUBROUTINE hsmt_spinor_soc

  
END MODULE m_hsmt_spinor
