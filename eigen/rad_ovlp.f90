!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

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
    REAL,    INTENT  (IN):: epar(0:,:,:)!(0:atoms%lmaxd,atoms%ntype,input%jspins)
    REAL,    INTENT  (IN):: vr(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
    REAL,    INTENT (OUT):: uun21(0:atoms%lmaxd,atoms%ntype),udn21(0:atoms%lmaxd,atoms%ntype)
    REAL,    INTENT (OUT):: dun21(0:atoms%lmaxd,atoms%ntype),ddn21(0:atoms%lmaxd,atoms%ntype)
    !     ..
    !     .. Local Scalars ..
    INTEGER itype,l,ispin,noded,nodeu,i
    REAL    wronk
    LOGICAL l_hia
    !     ..
    !     .. Local Arrays ..
    REAL :: f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    REAL :: g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)
    REAL :: vrTmp(atoms%jmtd)
    !     ..
    DO itype = 1,atoms%ntype
       DO l = 0,atoms%lmax(itype)
          l_hia=.FALSE.
          DO i = atoms%n_u+1, atoms%n_u+atoms%n_hia
             IF(atoms%lda_u(i)%atomType.EQ.itype.AND.atoms%lda_u(i)%l.EQ.l) THEN
                l_hia=.TRUE.
             ENDIF
          ENDDO
          DO ispin = 1,input%jspins

             !In the case of a spin-polarized calculation with Hubbard 1 we want to treat 
             !the correlated orbitals with a non-spin-polarized basis
             IF(l_hia.AND.input%jspins.EQ.2) THEN
                vrTmp = (vr(:,0,iType,1) + vr(:,0,iType,2))/2.0 
             ELSE
                vrTmp = vr(:,0,iType,ispin)
             ENDIF

             CALL radfun(l,iType,ispin,epar(l,iType,ispin),vrTmp,atoms,&
                 f(1,1,l,ispin),g(1,1,l,ispin),usdus,nodeu,noded,wronk) 

          ENDDO
          CALL int_21(f,g,atoms,itype,l,uun21,udn21,dun21,ddn21)
       ENDDO
    ENDDO

  END SUBROUTINE rad_ovlp
END

