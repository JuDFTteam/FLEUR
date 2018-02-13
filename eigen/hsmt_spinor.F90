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
  SUBROUTINE hsmt_spinor(isp,n,noco,chi_mat,chi_so)
    USE m_types
    IMPLICIT NONE

    TYPE(t_noco),INTENT(IN)      :: noco
    INTEGER,INTENT(IN)           :: isp, n
    COMPLEX,INTENT(OUT)          :: chi_mat(2,2)
    COMPLEX,INTENT(OUT),OPTIONAL :: chi_so(2,2)
   
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

    IF (.not.PRESENT(chi_so)) RETURN
    !In the first variation SOC case the off-diagonal spinors are needed
    IF (isp==1) THEN
       isp1=2;isp2=1
    ELSE
       isp1=1;isp2=2
    END IF

    chi_so(1,1) = chi(1,isp1)*CONJG(chi(1,isp2))
    chi_so(2,1) = chi(2,isp1)*CONJG(chi(1,isp2))
    chi_so(2,2) = chi(2,isp1)*CONJG(chi(2,isp2))
    chi_so(1,2) = chi(1,isp1)*CONJG(chi(2,isp2))

  END SUBROUTINE hsmt_spinor


END MODULE m_hsmt_spinor
