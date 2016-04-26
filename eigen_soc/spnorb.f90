MODULE m_spnorb
  !*********************************************************************
  !     calls soinit to calculate the radial spin-orbit matrix elements:
  !     rsopp,rsopdpd,rsoppd,rsopdp
  !     and sets up the so - angular matrix elements (soangl)
  !     using the functions anglso and sgml.
  !*********************************************************************
CONTAINS
  SUBROUTINE spnorb(atoms,noco,input,mpi, enpara, vr,spav, rsopp,rsoppd,rsopdp,rsopdpd,&
       usdus, rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop, soangl)

    USE m_anglso
    USE m_sgml
    USE m_soinit 
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)   :: mpi

    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_usdus),INTENT(OUT)   :: usdus
    !     ..
    !     .. Scalar Arguments ..
    LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    REAL,    INTENT (OUT) :: rsopp  (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsoppd (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsopdp (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsopdpd(atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsoplop (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsoplopd(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsopdplo(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsopplo (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsoploplop(atoms%ntypd,atoms%nlod,atoms%nlod,2,2)
    COMPLEX, INTENT (OUT) :: soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,&
         atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2)
    !     ..
    !     .. Local Scalars ..
    INTEGER is1,is2,jspin1,jspin2,l,l1,l2,m1,m2,n
    !     ..
    !     .. Local Arrays ..
    INTEGER ispjsp(2)
    !     ..
    !     ..
    DATA ispjsp/1,-1/

    CALL soinit(atoms,input,enpara, vr,spav, rsopp,rsoppd,rsopdp,rsopdpd,&
         usdus, rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop)
    !
    IF (mpi%irank.EQ.0) THEN
       DO n = 1,atoms%ntype
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,1),l=1,3)
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,1),l=1,3)
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,1),l=1,3)
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,1),l=1,3)
       ENDDO
    ENDIF
8000 FORMAT (' spin - orbit parameter HR  ')
8001 FORMAT (8f8.4)
9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')
    !

    IF ((ABS(noco%theta).LT.0.00001).AND.(ABS(noco%phi).LT.0.00001)) THEN
       !
       !       TEST for real function sgml(l1,m1,is1,l2,m2,is2)
       !
       DO l1 = 1,atoms%lmaxd
          DO l2 = 1,atoms%lmaxd
             DO jspin1 = 1,2
                DO jspin2 = 1,2
                   is1=ispjsp(jspin1)
                   is2=ispjsp(jspin2)
                   DO m1 = -l1,l1,1
                      DO m2 = -l2,l2,1
                         soangl(l1,m1,jspin1,l2,m2,jspin2) =&
                              CMPLX(sgml(l1,m1,is1,l2,m2,is2),0.0)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSE
       !
       !       TEST for complex function anglso(teta,phi,l1,m1,is1,l2,m2,is2)
       ! 
       DO l1 = 1,atoms%lmaxd
          DO l2 = 1,atoms%lmaxd
             DO jspin1 = 1,2
                DO jspin2 = 1,2
                   is1=ispjsp(jspin1)
                   is2=ispjsp(jspin2)
                   !
                   DO m1 = -l1,l1,1
                      DO m2 = -l2,l2,1
                         soangl(l1,m1,jspin1,l2,m2,jspin2) =&
                              anglso(noco%theta,noco%phi,l1,m1,is1,l2,m2,is2)
                      ENDDO
                   ENDDO
                   !
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !
    ENDIF

    IF (mpi%irank.EQ.0) THEN
       WRITE (6,FMT=8002)
       DO jspin1 = 1,2
          DO jspin2 = 1,2
             WRITE (6,FMT=*) 'd-states:is1=',jspin1,',is2=',jspin2
             WRITE (6,FMT='(7x,7i8)') (m1,m1=-3,3,1)
             WRITE (6,FMT=8003) (m2, (soangl(3,m1,jspin1,3,m2,jspin2),&
                  m1=-3,3,1),m2=-3,3,1)
          ENDDO
       ENDDO
    ENDIF

8002 FORMAT (' so - angular matrix elements ')
8003 FORMAT (i8,14f8.4)

  END SUBROUTINE spnorb
END MODULE m_spnorb
