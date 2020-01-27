!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_tetraef
      use m_juDFT
! -----------------------------------------------------------------------
! This subroutine evaluates the density of states with the tetrahedron method
! and sets the weight factors needed for the charge density for bulk systems.
! Adapted for the FLEUR code                                          GB 2000
! -----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE tetra_ef(
     >                    jspins,nkpt,
     >                    lb,ub,eig,zc,xfac,
     >                    ntetra,itetra,voltet,
     <                    efermi,w)
c
      IMPLICIT NONE
c
C     ..Scalar Arguments ..
      INTEGER, INTENT (IN)  :: jspins,nkpt,ntetra
      REAL,    INTENT (IN)  :: lb,ub,zc,xfac
      REAL,    INTENT (OUT) :: efermi
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)    :: itetra(:,:) !(4,6*nkptd)
      REAL,    INTENT (IN)    :: voltet(:) !(6*nkpt)
      REAL,    INTENT (OUT)   ::   w(:,:,:) !(neigd,nkptd,jspd)
      REAL,    INTENT (INOUT) :: eig(:,:,:) !(neigd,nkptd,jspd)
C     ..
C     .. Local Variables ..
      INTEGER i,j,jspin,neig,nk,nelec,ncr,ntet,it
      REAL    elow,dlow,eup,dup,ttt,dfermi,wgs
C     ..
C     .. Local Arrays ..
      REAL weight(4),ecmax(2,size(w,1))
      REAL wght(2,size(w,2),size(w,1)),eval(4)
C     ..


      DO neig = 1,size(w,1)
        DO jspin = 1,jspins
          ecmax(jspin,neig) = -1.0e25
          DO nk = 1,nkpt
            wght(jspin,nk,neig) = 0.0e0
               w(neig,nk,jspin) = 0.0e0
            IF ( eig(neig,nk,jspin).GT.ecmax(jspin,neig) )
     $              ecmax(jspin,neig) = eig(neig,nk,jspin)
          ENDDO
        ENDDO
      ENDDO
c
c  check for energy degeneracies in tetrahedrons
c
      DO jspin = 1,jspins
        DO ntet = 1,ntetra
          DO neig = 1,size(w,1)
            DO i = 1,3
              DO j = i+1,4
                IF (abs(eig(neig,itetra(i,ntet),jspin) -
     +                  eig(neig,itetra(j,ntet),jspin)).LT.1.0e-7) THEN
                  eig(neig,itetra(i,ntet),jspin) =
     +            eig(neig,itetra(i,ntet),jspin) + 1.0e-7
                  eig(neig,itetra(j,ntet),jspin) =
     +            eig(neig,itetra(j,ntet),jspin) - 1.0e-7
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c
c calculate weight factors
c
      DO ntet=1,ntetra
        DO neig=1,size(w,1)
          DO jspin=1,jspins
            DO i=1,4
              eval(i) = eig(neig,itetra(i,ntet),jspin)   
            ENDDO
            IF (max(eval(1),eval(2),eval(3),eval(4)).LT.9999.9) THEN
              DO i=1,4
                weight(i) = 1.0
                DO j=1,4
                  IF (i.NE.j) weight(i) = weight(i)*(eval(j)-eval(i))
                ENDDO
              ENDDO
              DO i=1,4
                weight(i) = 6.0*voltet(ntet)/weight(i)
                wght(jspin,itetra(i,ntet),neig) =
     +          wght(jspin,itetra(i,ntet),neig) + weight(i)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c
!      xfac = 2.0/jspins
      DO neig = 1,size(w,1)
        DO nk = 1,nkpt
          DO jspin = 1,jspins
            wght(jspin,nk,neig)=xfac*wght(jspin,nk,neig)
          ENDDO
        ENDDO
      ENDDO
c
c---------------------------------------------------
c determine fermi energy
c---------------------------------------------------
c
      nelec = zc                                     ! determine # of electrons
      ncr   = 0                                      ! (modified gb2000)

      elow = lb                                      ! determine lower bound 
      dlow = ncr
      DO nk = 1,nkpt
        DO neig = 1,size(w,1)
          DO jspin = 1,jspins
            ttt = elow - eig(neig,nk,jspin)
              IF ( elow.GT.ecmax(jspin,neig) )
     +          ttt = ecmax(jspin,neig) - eig(neig,nk,jspin)
              IF (ttt.LT.0.0e0) ttt = 0.0e0
              dlow = dlow + wght(jspin,nk,neig)*ttt*ttt*ttt/6
          ENDDO
        ENDDO
      ENDDO
      IF (dlow.GT.nelec) THEN    
        WRITE (6,180) elow,dlow,nelec
        CALL juDFT_error("dos: valence band too high ",calledby
     +       ="tetra_ef")
      ENDIF
  180 FORMAT (' valence band too high ',/,
     + '  elow ',f10.5,' dlow ',f10.5,' nelec ',i5)

      it  = 0
      eup = ub                                      ! determine upper bound
 10   dup = ncr
      DO nk = 1,nkpt
        DO neig = 1,size(w,1)
          DO jspin = 1,jspins
            ttt = eup - eig(neig,nk,jspin)  
            IF ( eup.GT.ecmax(jspin,neig) )
     +         ttt = ecmax(jspin,neig) - eig(neig,nk,jspin)  
            IF (ttt.LT.0.0e0) ttt = 0.0e0
            dup = dup + wght(jspin,nk,neig)*ttt*ttt*ttt/6
          ENDDO
        ENDDO
      ENDDO
      IF ( (dup-nelec).LT.0.00001 ) THEN 
        eup = eup + 0.2
        it  = it + 1
        IF( it .gt. 10 ) THEN
          WRITE (6,200) eup,dup,nelec
          CALL juDFT_error("dos: valence band too low ",
     +               calledby ="tetra_ef")
        END IF
        GOTO 10
      ENDIF
  200 FORMAT (' valence band too low ',/,
     + '  eup  ',f10.5,' dup  ',f10.5,' nelec ',i5)

      DO WHILE ( (eup-elow).GT.1.0e-10 )          ! iterate for fermi-energy
        efermi = 0.5*(elow+eup)
        dfermi = real(ncr)
        DO nk = 1,nkpt
          DO neig = 1,size(w,1)
            DO jspin = 1,jspins
              ttt  =efermi-eig(neig,nk,jspin)  
              IF ( efermi.GT.ecmax(jspin,neig) )
     +              ttt = ecmax(jspin,neig) - eig(neig,nk,jspin)  
              IF (ttt.LT.0.0e0) ttt = 0.0e0
              dfermi = dfermi + wght(jspin,nk,neig)*ttt*ttt*ttt/6
            ENDDO
          ENDDO
        ENDDO
        IF (dfermi.GT.nelec) THEN
          eup = efermi
        ELSE
          elow = efermi
        ENDIF
      ENDDO
      WRITE (6,220) efermi,dfermi,nelec
  220 FORMAT (//,'>>> D O S <<<',//,'   fermi energy ',f10.5,
     +       ' dtot ',f10.5,' nelec ',i5)
c
c---------------------------------------------------
c calculate weight factors for charge density integration
c---------------------------------------------------
c
      DO ntet = 1,ntetra
        DO neig = 1,size(w,1)
          DO jspin = 1,jspins
            DO i=1,4
              eval(i)=eig(neig,itetra(i,ntet),jspin)
            ENDDO
            IF (max(eval(1),eval(2),eval(3),eval(4)).LT.9999.9) THEN

              DO i = 1,4
                weight(i) = 1.0
                DO j = 1,4
                  IF (i.NE.j) THEN
                    weight(i) = weight(i) * (eval(j) - eval(i))
                  ENDIF
                ENDDO
                weight(i) = 6.0 * voltet(ntet) / weight(i)
              ENDDO

              wgs = 0.0e0
              DO i = 1,4
                ttt = efermi - eval(i)
                IF (efermi.GT.ecmax(jspin,neig)) 
     +                     ttt = ecmax(jspin,neig) - eval(i)
                IF ( ttt.LT.0.0e0 ) ttt = 0.0e0
                wgs = wgs + ttt**3*weight(i)
              ENDDO
              wgs = wgs / 24.0

              DO i = 1,4
                w(neig,itetra(i,ntet),jspin) =
     $          w(neig,itetra(i,ntet),jspin) + wgs
              ENDDO

            ENDIF
          ENDDO
        ENDDO
      ENDDO

c      DO jspin = 1,jspins
c        DO nk = 1,nkpt
c          DO neig = 1,size(w,1)
c            w(neig,nk,jspin) = xfac * w(neig,nk,jspin)
c          ENDDO
c        ENDDO
c      ENDDO

      RETURN
      END SUBROUTINE tetra_ef
      END MODULE m_tetraef
