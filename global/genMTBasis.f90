MODULE m_genMTBasis

CONTAINS

SUBROUTINE genMTBasis(atoms,enpara,vTot,mpi,iType,jspin,l_write,usdus,f,g,flo)

   USE m_types
   USE m_radfun
   USE m_radflo

   IMPLICIT NONE

   TYPE(t_atoms),  INTENT(IN)    :: atoms
   TYPE(t_enpara), INTENT(IN)    :: enpara
   TYPE(t_potden), INTENT(IN)    :: vTot
   TYPE(t_mpi),    INTENT(IN)    :: mpi
   TYPE(t_usdus),  INTENT(INOUT) :: usdus

   INTEGER,        INTENT(IN)    :: iType
   INTEGER,        INTENT(IN)    :: jspin
   LOGICAL,        INTENT(IN)    :: l_write

   REAL,           INTENT(INOUT) :: f(atoms%jmtd,2,0:atoms%lmaxd)
   REAL,           INTENT(INOUT) :: g(atoms%jmtd,2,0:atoms%lmaxd)
   REAL,           INTENT(INOUT) :: flo(atoms%jmtd,2,atoms%nlod)


   INTEGER                       :: l,nodeu,noded
   REAL                          :: wronk

   IF (l_write) WRITE (6,FMT=8000) iType

   DO l = 0,atoms%lmax(iType)
      CALL radfun(l,iType,jspin,enpara%el0(l,iType,jspin),vTot%mt(:,0,iType,jspin),atoms,&
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

8000  FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',&
                      /,t32,'radial function',t79,'energy derivative',/,t3,&
                      'l',t8,'energy',t26,'value',t39,'derivative',t53,&
                      'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,&
                      'norm',t119,'wronskian')
8010  FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)

END SUBROUTINE genMTBasis

END MODULE m_genMTBasis
