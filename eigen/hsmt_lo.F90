!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_lo
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_lo(input,atoms,sym,cell,mpi,noco,lapw,ud,tlmplm,fj,gj,n,chi,isp,iintsp,jintsp,hmat,smat)
    USE m_hlomat
    USE m_slomat
    USE m_setabc1lo
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: ud
    TYPE(t_tlmplm),INTENT(IN)   :: tlmplm
    
    TYPE(t_lapwmat),INTENT(INOUT)::hmat,smat
    
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN)   :: n
    INTEGER, INTENT (IN) :: isp,iintsp,jintsp !spins
    COMPLEX, INTENT(IN)  :: chi
    
    !Arrays
    REAL,INTENT(IN)      :: fj(:,:,:),gj(:,:,:)
    !     ..
    !     .. Local Scalars ..
    INTEGER na,nn
    !     ..
    !     .. Local Arrays ..
    REAL alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod)

    na = sum(atoms%neq(:n-1))
    DO nn = 1,atoms%neq(n)
       na = na + 1
       IF ((atoms%invsat(na).EQ.0) .OR. (atoms%invsat(na).EQ.1)) THEN
          
          
          IF (atoms%nlo(n).GE.1) THEN
             !--->          set up the a,b and c  coefficients
             !--->          for the local orbitals, if necessary.
             !--->          actually, these are the fj,gj equivalents
             CALL setabc1lo(atoms,n,ud,isp, alo1,blo1,clo1) 
             
             !--->       add the local orbital contribution to the overlap and
             !--->       hamiltonian matrix, if they are used for this atom.
             
             CALL slomat(&
                  input,atoms,mpi,lapw,cell,noco,n,na,&
                  isp,ud, alo1,blo1,clo1,fj,gj,&
                  iintsp,jintsp,chi,smat)
             CALL hlomat(input,atoms,mpi,lapw,ud,tlmplm,sym,cell,noco,isp,&
                  n,na,fj,gj,alo1,blo1,clo1,iintsp,jintsp,chi,hmat)
          ENDIF
    END IF
    !--->    end loop over equivalent atoms
 END DO
 RETURN
END SUBROUTINE hsmt_lo
END MODULE m_hsmt_lo
