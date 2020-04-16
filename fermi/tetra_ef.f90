!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_tetraef
   ! -----------------------------------------------------------------------
   ! This subroutine evaluates the density of states with the tetrahedron method
   ! and sets the weight factors needed for the charge density for bulk systems.
   ! Adapted for the FLEUR code                                          GB 2000
   ! -----------------------------------------------------------------------
   USE m_types
   USE m_constants
   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE tetra_ef(jspins,nkpt,lb,ub,eig,zc,xfac,ntetra,itetra,voltet,efermi,w)

      INTEGER, INTENT (IN)    :: jspins,nkpt,ntetra
      REAL,    INTENT (IN)    :: lb,ub,zc,xfac
      INTEGER, INTENT (IN)    :: itetra(:,:) !(4,6*nkptd)
      REAL,    INTENT (IN)    :: voltet(:)   !(6*nkpt)
      REAL,    INTENT (INOUT) :: eig(:,:,:)  !(neigd,nkptd,jspd)
      REAL,    INTENT (OUT)   :: w(:,:,:)    !(neigd,nkptd,jspd)
      REAL,    INTENT (OUT)   :: efermi

      INTEGER :: i,j,jspin,iBand,ikpt,nelec,ncr,itet,it,icorn,jcorn
      REAL    :: elow,dlow,eup,dup,ttt,dfermi,wgs
      REAL    :: weight(4),ecmax(2,size(w,1))
      REAL    :: wght(2,size(w,2),size(w,1)),eval(4)


      DO iBand = 1,size(w,1)
         DO jspin = 1,jspins
            ecmax(jspin,iBand) = -1.0e25
            DO ikpt = 1,nkpt
               wght(jspin,ikpt,iBand) = 0.0e0
               w(iBand,ikpt,jspin) = 0.0e0
               IF ( eig(iBand,ikpt,jspin).GT.ecmax(jspin,iBand)) ecmax(jspin,iBand) = eig(iBand,ikpt,jspin)
            ENDDO
         ENDDO
      ENDDO
      !
      !  check for energy degeneracies in tetrahedrons
      !
      DO jspin = 1,jspins
         DO itet = 1,ntetra
            DO iBand = 1,size(w,1)
               DO i = 1,3
                  icorn = itetra(i,itet)
                  DO j = i+1,4
                     jcorn = itetra(j,itet)
                     IF (abs(eig(iBand,icorn,jspin)-eig(iBand,jcorn,jspin)).LT.1.0e-7) THEN
                        eig(iBand,icorn,jspin) = eig(iBand,icorn,jspin) + 1.0e-7
                        eig(iBand,jcorn,jspin) = eig(iBand,jcorn,jspin) - 1.0e-7
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      ! calculate weight factors
      !
      DO itet=1,ntetra
         DO iBand=1,size(w,1)
            DO jspin=1,jspins

               eval = eig(iBand,itetra(:,itet),jspin)

               IF (ANY(eval.GE.9999.9)) CYCLE

               DO i=1,4
                  weight(i) = 1.0
                  DO j=1,4
                     IF (i.NE.j) weight(i) = weight(i)*(eval(j)-eval(i))
                  ENDDO
               ENDDO
               DO i=1,4
                  icorn = itetra(i,itet)
                  weight(i) = 6.0*voltet(itet)/weight(i)
                  wght(jspin,icorn,iBand) = wght(jspin,icorn,iBand) + weight(i)
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !
      !xfac = 2.0/jspins
      DO iBand = 1,size(w,1)
         DO ikpt = 1,nkpt
            DO jspin = 1,jspins
               wght(jspin,ikpt,iBand)=xfac*wght(jspin,ikpt,iBand)
            ENDDO
         ENDDO
      ENDDO
      !
      !---------------------------------------------------
      ! determine fermi energy
      !---------------------------------------------------
      !
      nelec = zc                                     ! determine # of electrons
      ncr   = 0                                      ! (modified gb2000)

      elow = lb                                      ! determine lower bound
      dlow = ncr
      DO ikpt = 1,nkpt
         DO iBand = 1,size(w,1)
            DO jspin = 1,jspins
               ttt = elow - eig(iBand,ikpt,jspin)
               IF ( elow.GT.ecmax(jspin,iBand) ) ttt = ecmax(jspin,iBand) - eig(iBand,ikpt,jspin)
               IF (ttt.LT.0.0e0)                ttt = 0.0e0
               dlow = dlow + wght(jspin,ikpt,iBand)*ttt*ttt*ttt/6
            ENDDO
         ENDDO
      ENDDO
      IF (dlow.GT.nelec) THEN
        WRITE (oUnit,180) elow,dlow,nelec
        CALL juDFT_error("dos: valence band too high ",calledby="tetra_ef")
      ENDIF
180   FORMAT (' valence band too high ',/,&
              '  elow ',f10.5,' dlow ',f10.5,' nelec ',i5)

      it  = 0
      eup = ub                                      ! determine upper bound
 10   dup = ncr
      DO ikpt = 1,nkpt
         DO iBand = 1,size(w,1)
            DO jspin = 1,jspins
               ttt = eup - eig(iBand,ikpt,jspin)
               IF ( eup.GT.ecmax(jspin,iBand) ) ttt = ecmax(jspin,iBand) - eig(iBand,ikpt,jspin)
               IF (ttt.LT.0.0e0)               ttt = 0.0e0
               dup = dup + wght(jspin,ikpt,iBand)*ttt*ttt*ttt/6
            ENDDO
         ENDDO
      ENDDO

      IF ( (dup-nelec).LT.0.00001 ) THEN
         eup = eup + 0.2
         it  = it + 1
         IF( it .gt. 10 ) THEN
            WRITE (oUnit,200) eup,dup,nelec
            CALL juDFT_error("dos: valence band too low ",calledby ="tetra_ef")
         END IF
         GOTO 10
      ENDIF

  200 FORMAT (' valence band too low ',/,&
              '  eup  ',f10.5,' dup  ',f10.5,' nelec ',i5)

      DO WHILE ( (eup-elow).GT.1.0e-10 )          ! iterate for fermi-energy
         efermi = 0.5*(elow+eup)
         dfermi = real(ncr)
         DO ikpt = 1,nkpt
            DO iBand = 1,size(w,1)
               DO jspin = 1,jspins
                  ttt  =efermi-eig(iBand,ikpt,jspin)
                  IF ( efermi.GT.ecmax(jspin,iBand) ) ttt = ecmax(jspin,iBand) - eig(iBand,ikpt,jspin)
                  IF (ttt.LT.0.0e0)                  ttt = 0.0e0
                  dfermi = dfermi + wght(jspin,ikpt,iBand)*ttt*ttt*ttt/6
               ENDDO
            ENDDO
         ENDDO
         IF (dfermi.GT.nelec) THEN
            eup = efermi
         ELSE
            elow = efermi
         ENDIF
      ENDDO

      WRITE (oUnit,220) efermi,dfermi,nelec
220   FORMAT (//,'>>> D O S <<<',//,'   fermi energy =',f10.5,&
             ' dtot =',f10.5,' nelec =',i5)
      !
      !---------------------------------------------------
      ! calculate weight factors for charge density integration
      !---------------------------------------------------
      !
      DO itet = 1,ntetra
         DO iBand = 1,size(w,1)
            DO jspin = 1,jspins
               eval = eig(iBand,itetra(:,itet),jspin)

               IF (ANY(eval.GE.9999.9)) CYCLE

               DO i = 1,4
                  weight(i) = 1.0
                  DO j = 1,4
                     IF (i.NE.j) THEN
                        weight(i) = weight(i) * (eval(j) - eval(i))
                     ENDIF
                  ENDDO
                  weight(i) = 6.0 * voltet(itet) / weight(i)
               ENDDO

               wgs = 0.0e0
               DO i = 1,4
                  ttt = efermi - eval(i)
                  IF (efermi.GT.ecmax(jspin,iBand)) ttt = ecmax(jspin,iBand) - eval(i)
                  IF ( ttt.LT.0.0e0 )              ttt = 0.0e0
                  wgs = wgs + ttt**3*weight(i)
               ENDDO
               wgs = wgs / 24.0

               w(iBand,itetra(:,itet),jspin) = w(iBand,itetra(:,itet),jspin) + wgs

            ENDDO
         ENDDO
      ENDDO

!     DO jspin = 1,jspins
!        DO ikpt = 1,nkpt
!           DO iBand = 1,size(w,1)
!              w(iBand,ikpt,jspin) = xfac * w(iBand,ikpt,jspin)
!           ENDDO
!        ENDDO
!     ENDDO

   END SUBROUTINE tetra_ef
END MODULE m_tetraef
