MODULE m_uj2f
  USE m_juDFT
!  *********************************************************************
!  * The calculation of slater integrals from u&j                      *
!  * input in eV; output in htr.                                       *
!  *-------------------------------------------------------------------*
!  * Extension to multiple U per atom type by G.M. 2017                *
!  *********************************************************************
CONTAINS
  SUBROUTINE uj2f(&
       jspins,atoms,&
       f0,f2,f4,f6)

    USE m_types
    IMPLICIT NONE
    !
    !  .. Arguments ..
    INTEGER,        INTENT(IN)  :: jspins
    TYPE(t_atoms),  INTENT(IN)  :: atoms
    REAL,           INTENT (OUT):: f0(atoms%n_u,jspins),f2(atoms%n_u,jspins)
    REAL,           INTENT (OUT):: f4(atoms%n_u,jspins),f6(atoms%n_u,jspins)
    !
    !  .. Local variables ..
    INTEGER l,itype,ltest,ispin,i_u
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
          DO i_u = 1, atoms%n_u
             itype = atoms%lda_u(i_u)%atomType
             l = atoms%lda_u(i_u)%l
             f2(i_u,ispin)=0.0 ; f4(i_u,ispin)=0.0 ; f6(i_u,ispin)=0.0
100          READ (45,'(i3,4f20.10)') ltest,ftest(1:4)
             IF (ltest.EQ.l) THEN
                f0(i_u,ispin) = ftest(1)
                IF (l.GT.0) THEN
                   f2(i_u,ispin) = ftest(2)
                   IF (l.GT.1) THEN
                      f4(i_u,ispin) = ftest(3)
                      IF (l.GT.2) THEN
                         f6(i_u,ispin) = ftest(4)
                      END IF
                   END IF
                END IF
             ELSE
                GOTO 100
             END IF
             READ (45,'(i3,4f20.10)') ltest,ftest(1)
             !                IF (ltest.EQ.0) THEN
             !                   f0(n,ispin) = f0(n,ispin) - ftest(1)
             !                ENDIF

             !              write(*,*) n,ispin,l,f0(n,ispin),f2(n,ispin),
             !    +                              f4(n,ispin),f6(n,ispin)
          END DO
       ENDDO
       CLOSE (45)
    ELSE
       !
       ! lda_u%l: orb.mom; lda_u%u,j: in eV
       !
       DO i_u = 1, atoms%n_u
          itype = atoms%lda_u(i_u)%atomType
          l = atoms%lda_u(i_u)%l
          u = atoms%lda_u(i_u)%u
          j = atoms%lda_u(i_u)%j
          !
          !        l.eq.0 :  f0 = u (the l=0 and l=1 case approximated g.b.`01)
          !
          IF (l.EQ.0) THEN
             f0(i_u,1) = u
             f2(i_u,1) = 0.0
             f4(i_u,1) = 0.0
             f6(i_u,1) = 0.0
             IF (j>0.00001) CALL juDFT_error("lda+u: no magnetic s-states", calledby ="uj2f")
             !
             !        l == 1 :  j = f2 / 5  (from PRL 80,5758 g.b.)
             !
          ELSE IF (l.EQ.1) THEN
             f0(i_u,1) = u
             f2(i_u,1) = 5.0*j
             f4(i_u,1) = 0.0
             f6(i_u,1) = 0.0
             !
             !        l.eq.2 : 3d: j=(f2+f4)/14; f4/f2 = 0.625
             !
          ELSE IF (l.EQ.2) THEN
             !             PRINT*, 'd-states'
             f0(i_u,1) = u
             f2(i_u,1) = 14.0*j/1.625
             f4(i_u,1) = f2(i_u,1)*0.625
             f6(i_u,1) = 0.0
             !
             !        l.eq. 3 : 4f: j=(286f2+195f4+250f6)/6435; f2/f4 = 675/451; f2/f6=2025/1001
             !
          ELSE IF (l.EQ.3) THEN
             !             PRINT*, 'f-states'
             f0(i_u,1) = u
             a= 286.0 + 195.0*451.0/675.0 + 250.0*1001.0/2025.0
             f2(i_u,1) = 6435.0*j/a
             f4(i_u,1) = 451.0/675.0*f2(i_u,1)
             f6(i_u,1) = 1001.0/2025.0*f2(i_u,1)
          ELSE
             PRINT*, 'lda+U is restricted to l<=3 ! You used l=', l
          END IF
          IF (jspins.EQ.2) THEN
             f0(i_u,jspins) = f0(i_u,1)
             f2(i_u,jspins) = f2(i_u,1)
             f4(i_u,jspins) = f4(i_u,1)
             f6(i_u,jspins) = f6(i_u,1)
          ENDIF

       END DO ! atoms%n_u
       ! 
    ENDIF

  END SUBROUTINE uj2f
END MODULE m_uj2f
