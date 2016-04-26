MODULE m_int21lo
  !-----------------------------------------------------------
  !
  ! Integrates combinations of flo(:,:,llo,s) with itself, f(:,:,l,s) 
  ! and g(:,:,l,s) (s=1,2) for given llo and itype with different 
  ! spins (s).
  ! Output is ..n21(l,itype), where .. is a combination of (u,d) and
  ! ulo dependent on the (f,g) combination used. Also ..n12 and 
  ! uloulopn21 are calculated.
  !
  !-----------------------------------------------------------
CONTAINS
  SUBROUTINE int_21lo(f,g,atoms,n, flo,ilo, uulon21,dulon21,uulon12,dulon12,uloulopn21)

    USE m_intgr, ONLY : intgr3
    USE m_types
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ilo,n
    !     ... Array Arguments
    REAL,    INTENT (IN) :: f(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,dimension%jspd)
    REAL,    INTENT (IN) :: g(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,dimension%jspd)
    REAL,    INTENT (IN) :: flo(:,:,:,:)!(atoms%jmtd,2,atoms%nlod,dimension%jspd)
    REAL,    INTENT (OUT):: uulon21,uulon12
    REAL,    INTENT (OUT):: dulon21,dulon12
    REAL,    INTENT (OUT):: uloulopn21(atoms%nlod,atoms%nlod)

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
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),uulon21)
    DO iri = 1, atoms%jri(n)
       uu_tmp(iri) = f(iri,1,l,1)*flo(iri,1,ilo,2)+ f(iri,2,l,1)*flo(iri,2,ilo,2)
    ENDDO
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),uulon12)
    !
    ! --> norm of product of du and ulo:
    !
    DO iri = 1, atoms%jri(n)
       uu_tmp(iri) = g(iri,1,l,2)*flo(iri,1,ilo,1) + g(iri,2,l,2)*flo(iri,2,ilo,1)
    ENDDO
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),dulon21)
    DO iri = 1, atoms%jri(n)
       uu_tmp(iri) = g(iri,1,l,1)*flo(iri,1,ilo,2) + g(iri,2,l,1)*flo(iri,2,ilo,2)
    ENDDO
    CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),dulon12)
    !
    ! --> norm of product of ulo and ulo':
    !
    DO ilop = 1, atoms%nlo(n)
       lp = atoms%llo(ilop,n)
       IF (l.EQ.lp) THEN
          DO iri = 1, atoms%jri(n)
             uu_tmp(iri) = flo(iri,1,ilo,2)*flo(iri,1,ilop,1) + flo(iri,2,ilo,2)*flo(iri,2,ilop,1)
          ENDDO
          CALL intgr3(uu_tmp,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),uloulopn21(ilo,ilop))

       ELSE
          uloulopn21(ilo,ilop) = 0.0
       ENDIF
    ENDDO

  END SUBROUTINE int_21lo
END MODULE m_int21lo
