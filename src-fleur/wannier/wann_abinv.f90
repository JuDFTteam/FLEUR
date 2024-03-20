!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_wann_abinv
CONTAINS
  SUBROUTINE wann_abinv(atoms,sym,acof,bcof,ccof)
    !     ***************************************************************
    !     Transform acof,bcof,ccof in case of atoms related by inversion
    !     symmetry to obtain the coefficients in the global frame.
    !     Based on abcrot.
    !     Frank Freimuth
    !     ***************************************************************
    USE m_types
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_sym),INTENT(IN)   :: sym

    COMPLEX, INTENT (INOUT) :: acof(:,0:,:)
    COMPLEX, INTENT (INOUT) :: bcof(:,0:,:)
    COMPLEX, INTENT (INOUT) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,neigd,nlod,natd)

    !     .. Local Scalars ..
    INTEGER :: itype,ineq,iatom,ilo,l

    iatom=0
    DO itype=1,atoms%ntype
       DO ineq=1,atoms%neq(itype)
          iatom=iatom+1
          IF(sym%invsat(iatom).NE.2) CYCLE
          DO l=1,atoms%lmax(itype),2
             acof(:,l**2:l*(l+2),iatom) = (-1)**l *&
                  acof(:,l**2:l*(l+2),iatom)
             bcof(:,l**2:l*(l+2),iatom) = (-1)**l * &
                  bcof(:,l**2:l*(l+2),iatom)
          ENDDO
          DO ilo=1,atoms%nlo(itype)
             l=atoms%llo(ilo,itype)
             IF(l.GT.0) THEN
                IF(MOD(l,2).EQ.0)CYCLE 
                ccof(-l:l,:,ilo,iatom) = (-1)**l * &
                     ccof(-l:l,:,ilo,iatom)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE wann_abinv
END MODULE m_wann_abinv
      

    
