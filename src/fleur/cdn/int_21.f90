MODULE m_int21
  !-----------------------------------------------------------
  !
  ! Integrates f(:,:,l,s) and g(:,:,l,s) (s=1,2) for given l
  ! and itype with different spins (s).
  ! Output is ..n21(l,itype), where .. is a (u,d) combination 
  ! dependet on the (f,g) combination used. 
  !
  !-----------------------------------------------------------
CONTAINS

  SUBROUTINE int_21(f,g,atoms,ityp,l,uu21n,ud21n,du21n,dd21n)
    
    USE m_intgr, ONLY : intgr3
    USE m_types_setup

    IMPLICIT NONE

    TYPE(t_atoms),            INTENT(IN)    :: atoms

    INTEGER, INTENT (IN) :: l,ityp

    REAL,    INTENT (IN) :: f(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    REAL,    INTENT (IN) :: g(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    REAL,    INTENT (INOUT) :: uu21n(0:atoms%lmaxd,atoms%ntype),ud21n(0:atoms%lmaxd,atoms%ntype)
    REAL,    INTENT (INOUT) :: du21n(0:atoms%lmaxd,atoms%ntype),dd21n(0:atoms%lmaxd,atoms%ntype)

    REAL        uu_tmp(atoms%jri(ityp))

    uu_tmp(:atoms%jri(ityp)) = f(:atoms%jri(ityp),1,l,2)*f(:atoms%jri(ityp),1,l,1)&
         + f(:atoms%jri(ityp),2,l,2)*f(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),uu21n(l,ityp))
    
    uu_tmp(:atoms%jri(ityp)) = f(:atoms%jri(ityp),1,l,2)*g(:atoms%jri(ityp),1,l,1)&
         + f(:atoms%jri(ityp),2,l,2)*g(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),ud21n(l,ityp))
    
    uu_tmp(:atoms%jri(ityp)) = g(:atoms%jri(ityp),1,l,2)*f(:atoms%jri(ityp),1,l,1)&
         + g(:atoms%jri(ityp),2,l,2)*f(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),du21n(l,ityp))
    
    uu_tmp(:atoms%jri(ityp)) = g(:atoms%jri(ityp),1,l,2)*g(:atoms%jri(ityp),1,l,1)&
         + g(:atoms%jri(ityp),2,l,2)*g(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),dd21n(l,ityp))
    
  END SUBROUTINE int_21

END MODULE m_int21
