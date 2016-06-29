!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_spinor
#include "juDFT_env.h"
  IMPLICIT NONE
CONTAINS

  !The spinors are calculated both in hssphn_sph & hssphn_nonsph, hence this is a
  !common subroutine
  SUBROUTINE hsmt_spinor(isp,n, noco,input,chi, chi11, chi21, chi22,l_socfirst,&
       isigma,isigma_x,isigma_y,isigma_z,chi11so,chi21so,chi22so,angso,chj,cross_k)
    USE m_types
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_noco),INTENT(IN)    :: noco
    INTEGER,INTENT(IN)  :: isp, n

    COMPLEX,INTENT(OUT) :: chi(2,2)
    COMPLEX,INTENT(OUT) :: chi11,chi21,chi22

    LOGICAL,INTENT(IN),OPTIONAL:: l_socfirst
    COMPLEX,INTENT(IN),OPTIONAL:: isigma(:,:,:)
    REAL,INTENT(IN),OPTIONAL   :: cross_k(:,:)
    COMPLEX,INTENT(OUT),OPTIONAL:: chj(:,:,:,:)
    COMPLEX,INTENT(OUT),OPTIONAL :: isigma_x(:,:),isigma_y(:,:),isigma_z(:,:)
    COMPLEX,INTENT(OUT),OPTIONAL :: chi11so(:,:),chi21so(:,:),chi22so(:,:)
    COMPLEX,INTENT(OUT),OPTIONAL :: angso(:,:,:)

    INTEGER           :: nsp,j1,j2,kj
    COMPLEX,PARAMETER :: ci=cmplx(0.0,1.0)

    !--->       set up the spinors of this atom within global
    !--->       spin-coordinateframe
    chi(1,1) =  exp(-ci*noco%alph(n)/2)*cos(noco%beta(n)/2)
    chi(1,2) = -exp(-ci*noco%alph(n)/2)*sin(noco%beta(n)/2)
    chi(2,1) =  exp(ci*noco%alph(n)/2)*sin(noco%beta(n)/2)
    chi(2,2) =  exp(ci*noco%alph(n)/2)*cos(noco%beta(n)/2)
    !--->       and determine the prefactors for the Hamitonian- and
    !--->       overlapp-matrix elements
    chi11 = chi(1,isp)*conjg(chi(1,isp))
    chi21 = chi(2,isp)*conjg(chi(1,isp))
    chi22 = chi(2,isp)*conjg(chi(2,isp))

    IF (.NOT.present(chj)) RETURN
    IF (noco%l_constr) THEN
       nsp = 3 - isp
       chj(isp,1,1,n) = chi(1,nsp)*conjg(chi(1,isp))
       chj(isp,2,1,n) = chi(2,nsp)*conjg(chi(1,isp))
       chj(isp,2,2,n) = chi(2,nsp)*conjg(chi(2,isp))
    ENDIF
    IF ((isp.EQ.2).AND.l_socfirst) THEN
       isigma_x=MATMUL(conjg(TRANSPOSE(chi)), MATMUL(isigma(:,:,1),chi))
       isigma_y=MATMUL(conjg(TRANSPOSE(chi)), MATMUL(isigma(:,:,2),chi))
       isigma_z=MATMUL(conjg(TRANSPOSE(chi)), MATMUL(isigma(:,:,3),chi))
       DO j1=1,input%jspins
          DO j2=1,input%jspins
             chi11so(j1,j2)=chi(1,j1)*conjg(chi(1,j2))
             chi21so(j1,j2)=chi(2,j1)*conjg(chi(1,j2))
             chi22so(j1,j2)=chi(2,j1)*conjg(chi(2,j2))
             DO kj = 1,size(angso,1)
                angso(kj,j1,j2)= isigma_x(j1,j2)*cross_k(kj,1)+&
                     isigma_y(j1,j2)*cross_k(kj,2)+ isigma_z(j1,j2)*cross_k(kj,3)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE hsmt_spinor


END MODULE m_hsmt_spinor
