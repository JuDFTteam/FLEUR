MODULE m_int21lo
  !-----------------------------------------------------------
  !
  ! Integrates combinations of flo(:,:,llo,s) with itself, f(:,:,l,s) 
  ! and g(:,:,l,s) (s=1,2) for given llo and itype with different 
  ! spins (s).
  ! Output is ..n21(l,itype), where .. is a combination of (u,d) and
  ! ulo dependent on the (f,g) combination used. Also ..n12 and 
  ! uloulop21n are calculated.
  !
  !-----------------------------------------------------------
CONTAINS
  SUBROUTINE int_21lo(f,g,atoms,n,flo,ilo,uulo21n,ulou21n,dulo21n,ulod21n,uloulop21n)

    USE m_intgr, ONLY : intgr3
    USE m_types_setup
    IMPLICIT NONE
    TYPE(t_atoms),             INTENT(IN)    :: atoms
    REAL,                      INTENT(INOUT) :: uulo21n(atoms%nlod,atoms%ntype)
    REAL,                      INTENT(INOUT) :: ulou21n(atoms%nlod,atoms%ntype)
    REAL,                      INTENT(INOUT) :: dulo21n(atoms%nlod,atoms%ntype)
    REAL,                      INTENT(INOUT) :: ulod21n(atoms%nlod,atoms%ntype)
    REAL,                      INTENT(INOUT) :: uloulop21n(atoms%nlod,atoms%nlod,atoms%ntype)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ilo,n
    !     ... Array Arguments
    REAL,    INTENT (IN) :: f(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    REAL,    INTENT (IN) :: g(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    REAL,    INTENT (IN) :: flo(:,:,:,:)!(atoms%jmtd,2,atoms%nlod,input%jspins)

    !     ...local scalars
    INTEGER iri,l,lp,ilop
    !     ...local arrays
    REAL    uu_tmp(atoms%jri(n))

    !
    ! --> norm of product of u and ulo:
    !
    l = atoms%llo(ilo,n)
    DO iri = 1, atoms%jri(n)
       uu_tmp(iri) = f(iri,1,l,2)*flo(iri,1,ilo,1)+ f(iri,2,l,2)*flo(iri,2,ilo,1)
    ENDDO
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),uulo21n(ilo,n))
    DO iri = 1, atoms%jri(n)
       uu_tmp(iri) = f(iri,1,l,1)*flo(iri,1,ilo,2)+ f(iri,2,l,1)*flo(iri,2,ilo,2)
    ENDDO
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),ulou21n(ilo,n))
    !
    ! --> norm of product of du and ulo:
    !
    DO iri = 1, atoms%jri(n)
       uu_tmp(iri) = g(iri,1,l,2)*flo(iri,1,ilo,1) + g(iri,2,l,2)*flo(iri,2,ilo,1)
    ENDDO
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),dulo21n(ilo,n))
    DO iri = 1, atoms%jri(n)
       uu_tmp(iri) = g(iri,1,l,1)*flo(iri,1,ilo,2) + g(iri,2,l,1)*flo(iri,2,ilo,2)
    ENDDO
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),ulod21n(ilo,n))
    !
    ! --> norm of product of ulo and ulo':
    !
    DO ilop = 1, atoms%nlo(n)
       lp = atoms%llo(ilop,n)
       IF (l.EQ.lp) THEN
          DO iri = 1, atoms%jri(n)
             uu_tmp(iri) = flo(iri,1,ilo,2)*flo(iri,1,ilop,1) + flo(iri,2,ilo,2)*flo(iri,2,ilop,1)
          ENDDO
          CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),uloulop21n(ilo,ilop,n))

       ELSE
          uloulop21n(ilo,ilop,n) = 0.0
       ENDIF
    ENDDO

  END SUBROUTINE int_21lo
END MODULE m_int21lo
