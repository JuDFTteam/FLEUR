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

   SUBROUTINE tetra_ef(kpts,jspins,lb,ub,eig,zc,xfac,efermi,w)

      TYPE(t_kpts),     INTENT(IN)    :: kpts
      INTEGER,          INTENT(IN)    :: jspins
      REAL,             INTENT(IN)    :: lb,ub,zc,xfac
      REAL,             INTENT(INOUT) :: eig(:,:,:)  !(neigd,nkptd,jspd)
      REAL,             INTENT(OUT)   :: w(:,:,:)    !(neigd,nkptd,jspd)
      REAL,             INTENT(OUT)   :: efermi

      INTEGER :: i,j,jspin,iBand,ikpt,nelec,ncr,itet,it,icorn,jcorn
      REAL    :: elow,dlow,eup,dup,ttt,dfermi,wgs
      REAL    :: weight(4),ecmax(2,size(w,1))
      REAL    :: wght(2,kpts%nkpt,size(w,1)),eval(4), tetra_weight(kpts%nkpt)

      w=0.0
      DO iBand = 1,size(w,1)
         DO jspin = 1,jspins
            ecmax(jspin,iBand) = -1.0e25
            DO ikpt = 1,kpts%nkpt
               wght(jspin,ikpt,iBand) = 0.0
               w(iBand,ikpt,jspin) = 0.0
               IF(eig(iBand,ikpt,jspin).GT.ecmax(jspin,iBand)) ecmax(jspin,iBand) = eig(iBand,ikpt,jspin)
            ENDDO
         ENDDO
      ENDDO
      !
      !  check for energy degeneracies in tetrahedrons
      !
      DO jspin = 1,jspins
         DO itet = 1,kpts%ntet
            DO iBand = 1,size(w,1)
               DO i = 1,3
                  icorn = kpts%ntetra(i,itet)
                  DO j = i+1,4
                     jcorn = kpts%ntetra(j,itet)
                     IF(abs(eig(iBand,icorn,jspin)-eig(iBand,jcorn,jspin)).LT.1.0e-7) THEN
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
      tetra_weight = 0.0
      DO itet=1,kpts%ntet
         DO i=1, 4
            tetra_weight(kpts%ntetra(i,itet)) = tetra_weight(kpts%ntetra(i,itet)) + kpts%voltet(itet)/(4.0*kpts%ntet)
         ENDDO
         DO iBand=1,size(w,1)
            DO jspin=1,jspins

               eval = eig(iBand,kpts%ntetra(:,itet),jspin)

               IF(ANY(eval.GE.9999.9)) CYCLE

               DO i=1,4
                  weight(i) = 1.0
                  DO j=1,4
                     IF(i.NE.j) weight(i) = weight(i)*(eval(j)-eval(i))
                  ENDDO
               ENDDO
               DO i=1,4
                  icorn = kpts%ntetra(i,itet)
                  weight(i) = 6.0*kpts%voltet(itet)/(weight(i)*kpts%ntet)
                  wght(jspin,icorn,iBand) = wght(jspin,icorn,iBand) + weight(i)
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !
      !xfac = 2.0/jspins
      tetra_weight = xfac*tetra_weight
      wght = xfac*wght
      !
      !---------------------------------------------------
      ! determine fermi energy
      !---------------------------------------------------
      !
      nelec = zc                                     ! determine # of electrons
      ncr   = 0                                      ! (modified gb2000)

      elow = lb                                      ! determine lower bound
      dlow = ncr
      DO ikpt = 1,kpts%nkpt
         DO iBand = 1,size(w,1)
            DO jspin = 1,jspins
               ttt = elow - eig(iBand,ikpt,jspin)
               IF ( elow.GT.ecmax(jspin,iBand) ) THEN
                  dlow = dlow + tetra_weight(ikpt)
                  cycle
               endif
               IF (ttt.LT.0.0e0) cycle
               dlow = dlow + wght(jspin,ikpt,iBand)*ttt*ttt*ttt/6
            ENDDO
         ENDDO
      ENDDO
      IF (dlow.GT.nelec) THEN
         WRITE (oUnit,180) elow,dlow,nelec
180      FORMAT (' valence band too high ',/,'  elow ',f10.5,' dlow ',f10.5,' nelec ',i5)
         CALL juDFT_error("dos: valence band too high ",calledby="tetra_ef")
      ENDIF


      it  = 0
      eup = ub
      dup = 0.0 ! determine upper bound
      DO WHILE ((dup-nelec).LT.0.00001)
         dup = ncr
         DO ikpt = 1,kpts%nkpt
            DO iBand = 1,size(w,1)
               DO jspin = 1,jspins
                  ttt = eup - eig(iBand,ikpt,jspin)
                  IF ( eup.GT.ecmax(jspin,iBand) ) THEN
                     dup = dup + tetra_weight(ikpt)
                     cycle
                  endif
                  IF (ttt.LT.0.0e0) cycle
                  dup = dup + wght(jspin,ikpt,iBand)*ttt*ttt*ttt/6
               ENDDO
            ENDDO
         ENDDO

         IF((dup-nelec).LT.0.00001) THEN
            eup = eup + 0.2
            it  = it + 1
            IF(it .gt. 10) THEN
               WRITE (oUnit,200) eup,dup,nelec
200            FORMAT (' valence band too low ',/,'  eup  ',f10.5,' dup  ',f10.5,' nelec ',i5)
               CALL juDFT_error("dos: valence band too low ",calledby ="tetra_ef")
            END IF
         ENDIF
      ENDDO



      DO WHILE ((eup-elow).GT.1.0e-10)          ! iterate for fermi-energy
         efermi = 0.5*(elow+eup)
         dfermi = real(ncr)
         DO ikpt = 1,kpts%nkpt
            DO iBand = 1,size(w,1)
               DO jspin = 1,jspins
                  ttt  =efermi-eig(iBand,ikpt,jspin)
                  IF ( efermi.GT.ecmax(jspin,iBand) ) THEN
                     dfermi = dfermi + tetra_weight(ikpt)
                     cycle
                  endif
                  IF (ttt.LT.0.0e0) cycle
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
220   FORMAT (//,'>>> D O S <<<',//,'   fermi energy =',f10.5,' dtot =',f10.5,' nelec =',i5)
      !
      !---------------------------------------------------
      ! calculate weight factors for charge density integration
      !---------------------------------------------------
      !
      DO itet = 1,kpts%ntet
         DO iBand = 1,size(w,1)
            DO jspin = 1,jspins
               eval = eig(iBand,kpts%ntetra(:,itet),jspin)

               IF(ANY(eval.GE.9999.9)) CYCLE

               DO i = 1,4
                  weight(i) = 1.0
                  DO j = 1,4
                     IF(i.NE.j) THEN
                        weight(i) = weight(i) * (eval(j) - eval(i))
                     ENDIF
                  ENDDO
                  weight(i) = 6.0 * kpts%voltet(itet)/(weight(i)*kpts%ntet)
               ENDDO

               IF ( efermi.GT.ecmax(jspin,iBand) ) THEN
                  wgs = kpts%voltet(itet)/(4.0* kpts%ntet)
               else
                  wgs = 0.0e0
                  DO i = 1,4
                     ttt = efermi - eval(i)
                     IF( ttt.LT.0.0e0 ) cycle
                     wgs = wgs + ttt**3*weight(i)
                  ENDDO
                  wgs = wgs / 24.0
               endif

               w(iBand,kpts%ntetra(:,itet),jspin) = w(iBand,kpts%ntetra(:,itet),jspin) + wgs

            ENDDO
         ENDDO
      ENDDO

!     DO jspin = 1,jspins
!        DO ikpt = 1,kpts%nkpt
!           DO iBand = 1,size(w,1)
!              w(iBand,ikpt,jspin) = xfac * w(iBand,ikpt,jspin)
!           ENDDO
!        ENDDO
!     ENDDO

   END SUBROUTINE tetra_ef
END MODULE m_tetraef
