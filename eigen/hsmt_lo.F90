!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_lo
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE Hsmt_lo(Input,Atoms,Sym,Cell,Mpi,Noco,nococonv,Lapw,Ud,Tlmplm,Fj,Gj,N,Chi,Isp,Iintsp,Jintsp,Hmat,Smat)
    USE m_hlomat
    USE m_slomat
    USE m_setabc1lo
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_nococonv),INTENT(IN) :: nococonv
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: ud
    TYPE(t_tlmplm),INTENT(IN)   :: tlmplm

    CLASS(t_mat),INTENT(INOUT)::hmat
    CLASS(t_mat),INTENT(INOUT),OPTIONAL::smat

    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN)   :: n
    INTEGER, INTENT (IN) :: isp,iintsp,jintsp !spins
    COMPLEX, INTENT(IN)  :: chi

    !Arrays
    REAL,INTENT(IN)      :: fj(:,:,:),gj(:,:,:)
    !     ..
    !     .. Local Scalars ..
    INTEGER na,nn,usp
    !     ..
    !     .. Local Arrays ..
    REAL alo1(atoms%nlod,input%jspins),blo1(atoms%nlod,input%jspins),clo1(atoms%nlod,input%jspins)
    CALL timestart("LO setup")

    na = SUM(atoms%neq(:n-1))
    DO nn = 1,atoms%neq(n)
       na = na + 1
       IF ((sym%invsat(na).EQ.0) .OR. (sym%invsat(na).EQ.1)) THEN


          IF (atoms%nlo(n).GE.1) THEN
             !--->          set up the a,b and c  coefficients
             !--->          for the local orbitals, if necessary.
             !--->          actually, these are the fj,gj equivalents
             DO usp=merge(1,isp,isp>2),merge(2,isp,isp>2)
               CALL setabc1lo(atoms,n,ud,usp,alo1,blo1,clo1)
             enddo

             !--->       add the local orbital contribution to the overlap and
             !--->       hamiltonian matrix, if they are used for this atom.

             IF (isp<3) THEN
                IF (.NOT.PRESENT(smat)) CALL judft_error("Bug in hsmt_lo, called without smat")
                CALL slomat(&
                     input,atoms,sym,mpi,lapw,cell,nococonv,n,na,&
                     isp,ud, alo1(:,isp),blo1(:,isp),clo1(:,isp),fj,gj,&
                     iintsp,jintsp,chi,smat)
             ENDIF
             CALL hlomat(input,atoms,mpi,lapw,ud,tlmplm,sym,cell,noco,nococonv,isp,&
                  n,na,fj,gj,alo1,blo1,clo1,iintsp,jintsp,chi,hmat)
          ENDIF
       END IF
       !--->    end loop over equivalent atoms
    END DO
    CALL timestop("LO setup")

    RETURN
  END SUBROUTINE hsmt_lo
END MODULE m_hsmt_lo
