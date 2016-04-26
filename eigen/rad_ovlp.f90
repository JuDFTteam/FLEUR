MODULE m_radovlp
  CONTAINS
  SUBROUTINE rad_ovlp(atoms,usdus,input,vr,epar, uun21,udn21,dun21,ddn21)
    !***********************************************************************
    ! calculates the overlap of the radial basis functions with different
    ! spin directions. These overlapp integrals are needed to calculate
    ! the contribution to the hamiltonian from the constant constraint
    ! B-field.
    !
    ! Philipp Kurz 2000-04
    !***********************************************************************

    USE m_int21
    USE m_radfun
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_usdus),INTENT(INOUT):: usdus

    !     .. Array Arguments ..
    REAL,    INTENT  (IN):: epar(0:,:,:)!(0:atoms%lmaxd,atoms%ntypd,dimension%jspd)
    REAL,    INTENT  (IN):: vr(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    REAL,    INTENT (OUT):: uun21(0:atoms%lmaxd,atoms%ntypd),udn21(0:atoms%lmaxd,atoms%ntypd)
    REAL,    INTENT (OUT):: dun21(0:atoms%lmaxd,atoms%ntypd),ddn21(0:atoms%lmaxd,atoms%ntypd)
    !     ..
    !     .. Local Scalars ..
    INTEGER itype,l,ispin,noded,nodeu
    REAL    wronk
    !     ..
    !     .. Local Arrays ..
    REAL :: f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    REAL :: g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    !     ..
    DO itype = 1,atoms%ntype
       DO l = 0,atoms%lmax(itype)
          DO ispin = 1,input%jspins
             CALL radfun(l,itype,ispin,epar(l,itype,ispin),vr(:,0,itype,ispin), atoms,&
                  f(1,1,l,ispin),g(1,1,l,ispin),usdus, nodeu,noded,wronk)
          ENDDO
          CALL int_21(f,g,atoms,itype,l, uun21(l,itype),udn21(l,itype), dun21(l,itype),ddn21(l,itype))
       ENDDO
    ENDDO

  END SUBROUTINE rad_ovlp
END

