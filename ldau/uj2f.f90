MODULE m_uj2f
  USE m_juDFT
!  *********************************************************************
!  * The calculation of slater integrals from u&j                      *
!  * input in eV; output in htr.                                       *
!  *********************************************************************
CONTAINS
  SUBROUTINE uj2f(&
       jspins,atoms,&
       f0,f2,f4,f6)

    USE m_types
    IMPLICIT NONE
    !
    !  .. Arguments ..
    INTEGER,INTENT(IN)       :: jspins
    TYPE(t_atoms),INTENT(IN) :: atoms
    REAL,           INTENT (OUT):: f0(atoms%n_u,jspins),f2(atoms%n_u,jspins)
    REAL,           INTENT (OUT):: f4(atoms%n_u,jspins),f6(atoms%n_u,jspins)
    !
    !  .. Local variables ..
    INTEGER n,l,itype,ltest,ispin
    REAL u,j,a,ftest(4)
    LOGICAL l_exist

    l_exist=.FALSE.
    INQUIRE (file='slaterf',exist=l_exist)

    IF (l_exist) THEN
       !
       ! --> f's have been calculated in cored ; read from file
       !
       OPEN (45,file='slaterf',form='formatted',status='old')
       DO ispin = 1, jspins
          n = 0
          DO itype=1,atoms%ntype
             l = atoms%lda_u(itype)%l
             IF (l.GE.0) THEN
                n = n + 1
                f2(n,ispin)=0.0 ; f4(n,ispin)=0.0 ; f6(n,ispin)=0.0
100             READ (45,'(i3,4f20.10)') ltest,ftest(1:4)
                IF (ltest.EQ.l) THEN
                   f0(n,ispin) = ftest(1)
                   IF (l.GT.0) THEN
                      f2(n,ispin) = ftest(2)
                      IF (l.GT.1) THEN
                         f4(n,ispin) = ftest(3)
                         IF (l.GT.2) THEN
                            f6(n,ispin) = ftest(4)
                         ENDIF
                      ENDIF
                   ENDIF
                ELSE
                   GOTO 100
                ENDIF
                READ (45,'(i3,4f20.10)') ltest,ftest(1)
                !                IF (ltest.EQ.0) THEN
                !                   f0(n,ispin) = f0(n,ispin) - ftest(1)
                !                ENDIF
             ENDIF
             !              write(*,*) n,ispin,l,f0(n,ispin),f2(n,ispin),
             !    +                              f4(n,ispin),f6(n,ispin)
          ENDDO
       ENDDO
       CLOSE (45)
    ELSE
       !
       ! lda_u%l: orb.mom; lda_u%u,j: in eV
       !
       n = 0
       DO itype=1,atoms%ntype
          l = atoms%lda_u(itype)%l
          IF (l.GE.0) THEN
             n = n + 1
             u = atoms%lda_u(itype)%u
             j = atoms%lda_u(itype)%j
             !
             !        l.eq.0 :  f0 = u (the l=0 and l=1 case approximated g.b.`01)
             !
             IF (l.EQ.0) THEN
                f0(n,1) = u
                f2(n,1) = 0.0
                f4(n,1) = 0.0
                f6(n,1) = 0.0
                IF (j>0.00001)  CALL juDFT_error&
                     &            ("lda+u: no magnetic s-states",calledby ="uj2f")
                !
                !        l == 1 :  j = f2 / 5  (from PRL 80,5758 g.b.)
                !
             ELSEIF (l.EQ.1) THEN
                f0(n,1) = u
                f2(n,1) = 5.0*j
                f4(n,1) = 0.0
                f6(n,1) = 0.0
                !
                !        l.eq.2 : 3d: j=(f2+f4)/14; f4/f2 = 0.625
                !
             ELSEIF (l.EQ.2) THEN
                !             PRINT*, 'd-states'
                f0(n,1) = u
                f2(n,1) = 14.0*j/1.625
                f4(n,1) = f2(n,1)*0.625
                f6(n,1) = 0.0
                !
                !        l.eq. 3 : 4f: j=(286f2+195f4+250f6)/6435; f2/f4 = 675/451; f2/f6=2025/1001
                !
             ELSEIF (l.EQ.3) THEN
                !             PRINT*, 'f-states'
                f0(n,1) = u
                a= 286.0 + 195.0*451.0/675.0 + 250.0*1001.0/2025.0
                f2(n,1) = 6435.0*j/a
                f4(n,1) = 451.0/675.0*f2(n,1)
                f6(n,1) = 1001.0/2025.0*f2(n,1)
             ELSE
                PRINT*, 'lda+U is restricted to l<=3 ! You used l=', l
             ENDIF
             IF (jspins.EQ.2) THEN
                f0(n,jspins) = f0(n,1)
                f2(n,jspins) = f2(n,1)
                f4(n,jspins) = f4(n,1)
                f6(n,jspins) = f6(n,1)
             ENDIF
          ENDIF
       ENDDO ! atoms%ntype
       ! 
    ENDIF

  END SUBROUTINE uj2f
END MODULE m_uj2f
