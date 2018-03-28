!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_soinit
  !
  !**********************************************************************
  !     generates radial spin-orbit matrix elements:sorad
  !**********************************************************************
  !
CONTAINS
  SUBROUTINE soinit(atoms,input,enpara, vr,spav,rsoc,usdus)

    USE m_sorad
    USE m_types
    IMPLICIT NONE
    TYPE(t_enpara),INTENT(IN)  :: enpara
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_usdus),INTENT(INOUT):: usdus
    TYPE(t_rsoc),INTENT(INOUT) :: rsoc
    !
    !     .. Scalar Arguments ..
    !     ..
    !     .. Scalar Arguments ..
    LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,jspin,n
    !     ..
    !     .. Local Arrays ..
    REAL  vr0(atoms%jmtd,size(vr,4))
    !     ..
    rsoc%rsopp  =0.0
    rsoppd =0.0
    rsopdp =0.0
    rsopdpd=0.0
    rsoplop =0.0
    rsoplopd=0.0
    rsopdplo=0.0
    rsopplo=0.0
    rsoploplop=0.0
    DO n = 1,atoms%ntype
       vr0=0.0
       vr0(:atoms%jri(n),:) = vr(:atoms%jri(n),0,n,:)
       
       CALL sorad(&
            atoms,input,n,vr0,enpara,spav,&
            rsopp,rsopdpd,rsoppd,rsopdp,usdus,&
            rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop)
       
    END DO ! end-do-loop : atoms%ntype

  END SUBROUTINE soinit
END MODULE m_soinit
