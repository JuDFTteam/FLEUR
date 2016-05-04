MODULE m_soinit
  !
  !**********************************************************************
  !     1. generates radial spin-orbit matrix elements:sorad
  !     2. generates spin-angular spin-orbit matrix   :soorb (not implemented)
  !**********************************************************************
  !
CONTAINS
  SUBROUTINE soinit(atoms,input,enpara, vr,spav,&
       rsopp,rsoppd,rsopdp,rsopdpd,usdus,&
       rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop)

    USE m_sorad
    USE m_types
    IMPLICIT NONE
    TYPE(t_enpara),INTENT(IN)  :: enpara
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_usdus),INTENT(INOUT)   :: usdus
    !
    !     .. Scalar Arguments ..
    !     ..
    !     .. Scalar Arguments ..
    LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    REAL,    INTENT (OUT) :: rsopp  (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsoppd (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsopdp (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsopdpd(atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsoplop (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsoplopd(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsopdplo(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsopplo (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsoploplop(atoms%ntypd,atoms%nlod,atoms%nlod,2,2)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,jspin,n
    !     ..
    !     .. Local Arrays ..
    REAL  vr0(atoms%jmtd,size(vr,4))
    !     ..
    rsopp  =0.0
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
