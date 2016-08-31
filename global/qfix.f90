      MODULE m_qfix
      USE m_juDFT
!     *******************************************************
!     check total charge and renormalize        c,l.fu
!     *******************************************************
      CONTAINS
      SUBROUTINE qfix(&
     &                stars,atoms,sym,vacuum,&
     &                sphhar,input,cell,oneD,&
     &                qpw,rhtxy,rho,rht,l_printData,&
     &                fix)

      USE m_types
      USE m_cdntot
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_stars),INTENT(IN) :: stars
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_vacuum),INTENT(IN):: vacuum
      TYPE(t_sphhar),INTENT(IN):: sphhar
      TYPE(t_input),INTENT(IN) :: input
      TYPE(t_oneD),INTENT(IN)  :: oneD
      TYPE(t_cell),INTENT(IN)  :: cell
      LOGICAL,INTENT(IN)       :: l_printData
      REAL,    INTENT (OUT) :: fix
!-odim
!+odim
!     ..
!     .. Array Arguments ..
      COMPLEX,INTENT (INOUT) :: qpw(stars%n3d,input%jspins)
      COMPLEX,INTENT (INOUT) :: rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,input%jspins)
      REAL,   INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins)
      REAL,   INTENT (INOUT) :: rht(vacuum%nmzd,2,input%jspins)
!     ..
!     .. Local Scalars ..
      LOGICAL fixtot,l99
      REAL qtot,qis,zc
      INTEGER ivac,j,jm,jspin,k,lh,n,nl,ns,na
!     ..
      CALL cdntot(&
     &            stars,atoms,sym,&
     &            vacuum,input,cell,oneD,&
     &            qpw,rho,rht,.FALSE.,&
     &            qtot,qis)
      zc = 0.
      DO 10 n = 1,atoms%ntype
         zc = zc + atoms%neq(n)*atoms%zatom(n)
   10 CONTINUE

!+roa (check via qfix file if total charge or only interstitial to fix)
      fixtot=.TRUE.
      INQUIRE(file='qfix',exist=l99)
      IF (l99) then
        OPEN (99,file='qfix',status='old',form='formatted')
        READ (99,'(1x,l1)',end=1199) fixtot
        IF (.NOT.fixtot ) THEN
           REWIND (99)
           WRITE (99,'(1x,l1,70a)') .TRUE.,&
     &      ' (1x,l1) F..fix interstitial T..fix total charge '
        ENDIF
 1199   CLOSE (99)
      ENDIF
      zc = zc + 2*input%efield%sigma
      IF ( fixtot ) THEN
!-roa
         fix = zc/qtot
         DO 100 jspin = 1,input%jspins
            na = 1
            DO 40 n = 1,atoms%ntype
               ns = atoms%ntypsy(na)
               lh = sphhar%nlh(ns)
               jm = atoms%jri(n)
               DO 30 nl = 0,lh
                  DO 20 j = 1,jm
                     rho(j,nl,n,jspin) = fix*rho(j,nl,n,jspin)
   20             CONTINUE
   30          CONTINUE
               na = na + atoms%neq(n)
   40       CONTINUE
            DO 50 k = 1,stars%ng3
               qpw(k,jspin) = fix*qpw(k,jspin)
   50       CONTINUE
            IF (input%film) THEN
               DO 90 ivac = 1,vacuum%nvac
                  DO 60 n = 1,vacuum%nmz
                     rht(n,ivac,jspin) = fix*rht(n,ivac,jspin)
   60             CONTINUE
                  DO 80 n = 1,vacuum%nmzxy
                     DO 70 k = 2,oneD%odi%nq2
                        rhtxy(n,k-1,ivac,jspin) = fix*&
     &                    rhtxy(n,k-1,ivac,jspin)
   70                CONTINUE
   80             CONTINUE
   90          CONTINUE
            END IF
  100    CONTINUE
         WRITE (6,FMT=8000) zc,fix
         CALL cdntot(&
     &               stars,atoms,sym,&
     &               vacuum,input,cell,oneD,&
     &               qpw,rho,rht,l_printData,&
     &               qtot,qis)
!+roa 
      ELSE
         fix = (zc - qtot) / qis + 1.
         DO jspin = 1,input%jspins
            DO k = 1,stars%ng3
               qpw(k,jspin) = fix*qpw(k,jspin)
            ENDDO
         ENDDO
         WRITE (6,FMT=8001) zc,fix
         CALL cdntot(&
     &               stars,atoms,sym,&
     &               vacuum,input,cell,oneD,&
     &               qpw,rho,rht,l_printData,&
     &               qtot,qis)

      ENDIF

      IF (fix>1.1) CALL juDFT_WARN("You lost too much charge")
      IF (fix<.9) CALL juDFT_WARN("You gained too much charge")


 8000 FORMAT (/,10x,'zc= ',f12.6,5x,'qfix=  ',f10.6)
 8001 FORMAT (/,' > broy only qis: ','zc= ',f12.6,5x,'qfix=  ',f10.6)
!-roa

      END SUBROUTINE qfix
      END MODULE m_qfix
