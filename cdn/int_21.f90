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
  SUBROUTINE int_21(f,g,atoms,ityp,l, uun21,udn21,dun21,ddn21)
    
    USE m_intgr, ONLY : intgr3
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: l,ityp
    REAL,    INTENT (OUT):: uun21,udn21,dun21,ddn21
    !     ... Array Arguments
    REAL,    INTENT (IN) :: f(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,dimension%jspd)
    REAL,    INTENT (IN) :: g(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,dimension%jspd)
    !     ...local arrays
    REAL        uu_tmp(atoms%jri(ityp))


    uu_tmp(:atoms%jri(ityp)) = f(:atoms%jri(ityp),1,l,2)*f(:atoms%jri(ityp),1,l,1)&
         + f(:atoms%jri(ityp),2,l,2)*f(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),uun21)
    
    uu_tmp(:atoms%jri(ityp)) = f(:atoms%jri(ityp),1,l,2)*g(:atoms%jri(ityp),1,l,1)&
         + f(:atoms%jri(ityp),2,l,2)*g(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),udn21)
    
    uu_tmp(:atoms%jri(ityp)) = g(:atoms%jri(ityp),1,l,2)*f(:atoms%jri(ityp),1,l,1)&
         + g(:atoms%jri(ityp),2,l,2)*f(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),dun21)
    
    uu_tmp(:atoms%jri(ityp)) = g(:atoms%jri(ityp),1,l,2)*g(:atoms%jri(ityp),1,l,1)&
         + g(:atoms%jri(ityp),2,l,2)*g(:atoms%jri(ityp),2,l,1)
    CALL intgr3(uu_tmp,atoms%rmsh(:,ityp),atoms%dx(ityp),atoms%jri(ityp),ddn21)
    
  END SUBROUTINE int_21
END MODULE m_int21
