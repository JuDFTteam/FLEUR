!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_genMTBasis

CONTAINS

  SUBROUTINE genMTBasis(atoms,enpara,vTot,mpi,iType,jspin,usdus,f,g,flo,l_dftspinpol)
    USE m_types
    USE m_radfun
    USE m_radflo
    !$  use omp_lib

    IMPLICIT NONE

    TYPE(t_atoms),  INTENT(IN)    :: atoms
    TYPE(t_enpara), INTENT(IN)    :: enpara
    TYPE(t_potden), INTENT(IN)    :: vTot
    TYPE(t_mpi),    INTENT(IN)    :: mpi
    TYPE(t_usdus),  INTENT(INOUT) :: usdus

    INTEGER,        INTENT(IN)    :: iType
    INTEGER,        INTENT(IN)    :: jspin

    LOGICAL,        INTENT(IN)    :: l_dftspinpol

    REAL,           INTENT(INOUT) :: f(atoms%jmtd,2,0:atoms%lmaxd)
    REAL,           INTENT(INOUT) :: g(atoms%jmtd,2,0:atoms%lmaxd)
    REAL,           INTENT(INOUT) :: flo(atoms%jmtd,2,atoms%nlod)


    INTEGER                       :: l,nodeu,noded
    REAL                          :: wronk


    LOGICAL    :: l_write,l_hia
    REAL       :: vrTmp(atoms%jmtd)
    INTEGER    :: i
    l_write=mpi%irank==0 
    !$ l_write = l_write .and. omp_get_num_threads()==1


    IF (l_write) WRITE (6,FMT=8000) iType

    DO l = 0,atoms%lmax(iType)
       !Check if the orbital is to be treated with Hubbard 1
       l_hia=.FALSE.
       DO i = atoms%n_u+1, atoms%n_u+atoms%n_hia
          IF(atoms%lda_u(i)%atomType.EQ.itype.AND.atoms%lda_u(i)%l.EQ.l) THEN
             l_hia=.TRUE.
          ENDIF
       ENDDO

       !In the case of a spin-polarized calculation with Hubbard 1 we want to treat 
       !the correlated orbitals with a non-spin-polarized basis  
       IF(l_hia.AND.SIZE(vTot%mt,4).GT.1.AND..NOT.l_dftspinpol) THEN
          vrTmp = (vTot%mt(:,0,iType,1) + vTot%mt(:,0,iType,2))/2.0
       ELSE
          vrTmp = vTot%mt(:,0,iType,jspin)
       ENDIF
       CALL radfun(l,iType,jspin,enpara%el0(l,iType,jspin),vrTmp,atoms,&
            f(1,1,l),g(1,1,l),usdus,nodeu,noded,wronk)      
       IF (l_write) THEN
          WRITE (6,FMT=8010) l,enpara%el0(l,iType,jspin),usdus%us(l,iType,jspin),usdus%dus(l,iType,jspin),&
               nodeu,usdus%uds(l,iType,jspin),usdus%duds(l,iType,jspin),noded,usdus%ddn(l,iType,jspin),wronk
       END IF
    END DO

    ! Generate the extra wavefunctions for the local orbitals, if there are any.
    IF (atoms%nlo(iType).GE.1) THEN
       CALL radflo(atoms,iType,jspin,enpara%ello0(1,1,jspin),vTot%mt(:,0,iType,jspin),f,g,mpi,&
            usdus,usdus%uuilon(1,1,jspin),usdus%duilon(1,1,jspin),usdus%ulouilopn(1,1,1,jspin),flo)
    END IF

8000 FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',&
         /,t32,'radial function',t79,'energy derivative',/,t3,&
         'l',t8,'energy',t26,'value',t39,'derivative',t53,&
         'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,&
         'norm',t119,'wronskian')
8010 FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)

  END SUBROUTINE genMTBasis

END MODULE m_genMTBasis
