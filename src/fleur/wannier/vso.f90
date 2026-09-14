!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vsoc
   implicit none
CONTAINS
  SUBROUTINE vsoc(input,atoms,vr,epar, l_spav, vso)
    !*************************************************
    !     Compute the spin-orbit potential.
    !*************************************************
    USE m_sointg
    USE m_types
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    REAL, INTENT(in)    :: vr(:,:,:)
    REAL, INTENT(in)    :: epar(:,:,:)
    LOGICAL,INTENT(in)  :: l_spav
    REAL, INTENT(out)   :: vso(:,:,:)

    REAL    :: v0(atoms%jmtd)
    REAL    :: e
    INTEGER :: n,i,l

    l=2

    DO n=1,atoms%ntype
       v0(:) = 0.0
       IF (input%jspins.EQ.1) THEN
          v0(1:atoms%jri(n)) = vr(1:atoms%jri(n),n,1)
          e = epar(l,n,1)
       ELSE
          DO i = 1,atoms%jri(n)
             v0(i) = (vr(i,n,1)+vr(i,n,input%jspins))/2.
          END DO
          e = (epar(l,n,1)+epar(l,n,input%jspins))/2.
       END IF

       CALL sointg(n,e,vr(:,n,:),v0,atoms,input, vso(:,n,:))  

       IF(l_spav)THEN
          DO i= 1,atoms%jmtd
             vso(i,n,1)= (vso(i,n,1)+vso(i,n,2))/2.
             vso(i,n,2)= vso(i,n,1)
          ENDDO
       ENDIF

    ENDDO



  END SUBROUTINE vsoc
END MODULE m_vsoc
