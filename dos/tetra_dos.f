!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_tetrados
c----------------------------------------------------------------------
c
c This subroutine evaluates the density of states (g) by the linear
c tetrahedron method on a energy grid (e) of 'ned' points.
c
c ev()          ... eigenvalues 
c qal()         ... partial charges
c ntetra)       ... number of tetrahedrons
c itetra(1-4,nt)... index of k-points forming tetrahedron nt
c voltet(nt)    ... volume of tetrahedron nt
c omega_bz      ... volume of irreducible part of BZ
c
c                                                      gb 2000
c----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE tetra_dos(
     >                     lmax,ntype,neigd,ned,ntetra,nkpt,
     >                     itetra,efermi,voltet,energy,nevk,
     >                     ev,qal,
     <                     g)

      IMPLICIT NONE
c
C     ..Scalar Arguments ..
      INTEGER, INTENT (IN)  :: ntype,neigd,ned,lmax
      INTEGER, INTENT (IN)  :: ntetra,nkpt
      REAL,    INTENT (IN)  :: efermi
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)    :: itetra(4,6*nkpt),nevk(nkpt)
      REAL,    INTENT (IN)    :: voltet(6*nkpt),energy(ned)
      REAL,    INTENT (IN)    :: qal(lmax*ntype+3,neigd,nkpt)
      REAL,    INTENT (INOUT) :: ev(neigd,nkpt)
      REAL,    INTENT (OUT)   :: g(ned,lmax*ntype+3)
C     ..
C     .. Local Variables ..
      INTEGER i,j,neig,nk,ntp,ne,ns,nc,ntet
      REAL    ener,efer,w
C     ..
C     .. Local Arrays ..
      REAL weight(4),eval(4),ecmax(neigd),term(ned)
      REAL wpar(4,ntype,neigd,nkpt),wparint(neigd,nkpt)
      REAL wght(neigd,nkpt)

      DO nk = 1,nkpt
         ev(nevk(nk)+1:neigd,nk) = 1.0e10 
      ENDDO
      wpar(:,:,:,:) = 0.0

      DO neig = 1,neigd
         ecmax(neig) = -1.0e25
         DO nk = 1,nkpt
            wght(neig,nk) = 0.0e0
            IF ( ev(neig,nk).GT.ecmax(neig) ) ecmax(neig) = ev(neig,nk)
         ENDDO
      ENDDO
c
c  check for energy degeneracies in tetrahedrons
c
      DO ntet = 1,ntetra
        DO neig = 1,neigd
          DO i = 1,3
            DO j = i+1,4
              IF (abs(ev(neig,itetra(i,ntet)) -
     +                ev(neig,itetra(j,ntet))).LT.1.0e-7) THEN
                   ev(neig,itetra(i,ntet))=
     +             ev(neig,itetra(i,ntet))+i*(1.e-7)*ntet
                   ev(neig,itetra(j,ntet))=
     +             ev(neig,itetra(j,ntet))-i*(1.e-7)*ntet
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE (*,*) 'in tetra_dos'  ! do not remove  this statement !
c
c calculate weight factors
c
      DO ntet = 1,ntetra
        w = 6.0*voltet(ntet)
        DO neig = 1,neigd
          DO i = 1,4
            eval(i) = ev(neig,itetra(i,ntet))
          ENDDO

          IF (max(eval(1),eval(2),eval(3),eval(4)).LT.9999.9) THEN
            DO i=1,4
              weight(i) = 1.0
              DO j=1,4
                IF (i.NE.j) weight(i) = weight(i)*(eval(j)-eval(i))
              ENDDO
              wght(neig,itetra(i,ntet)) =
     +        wght(neig,itetra(i,ntet)) + w/weight(i)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c
c calculate partial weights
c
      DO nk=1,nkpt
        DO neig = 1,nevk(nk)
          DO ntp = 1,ntype
            nc = lmax*(ntp-1)

            DO ntet = 1,ntetra
              IF (ALL(itetra(:,ntet).ne.nk)) CYCLE
           
              eval(1:4) = ev(neig,itetra(1:4,ntet))

              IF (max(eval(1),eval(2),eval(3),eval(4)).LT.9999.9) THEN
                DO i=1,4
                  weight(i)=1.0
                  DO j=1,4
                    IF (i.NE.j) weight(i)=weight(i)*(eval(j)-eval(i))
                  ENDDO
                  weight(i)=6.0*voltet(ntet)/weight(i)
                  DO ns=1,4
                    wpar(ns,ntp,neig,itetra(i,ntet)) =
     +              wpar(ns,ntp,neig,itetra(i,ntet)) + 
     +               0.25*weight(i)*qal(nc+ns,neig,nk)
                  ENDDO
                  IF (ntp.EQ.1) wparint(neig,itetra(i,ntet)) = 
     +                          wparint(neig,itetra(i,ntet)) + 
     +                0.25*weight(i)*qal(lmax*ntype+1,neig,nk)
                ENDDO
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO
c
c---------------------------------------------------
c evaluate d.o.s.
c---------------------------------------------------
c
      g(:,:) = 0.0

      DO nk = 1,nkpt
        DO neig = 1,neigd
          
          ener = ev(neig,nk)
          DO ntp = 1,ntype
            DO ns = 1,lmax
              nc = ns + lmax*(ntp-1)
              w  = 0.5*wpar(ns,ntp,neig,nk)
              DO ne = 1,ned
                term(ne) = energy(ne) - ener
                IF ( energy(ne).GT.ecmax(neig) )
     +                term(ne) = ecmax(neig) - ener
                IF ( term(ne).LT.0.0e0 ) term(ne) = 0.0e0
                g(ne,nc) = g(ne,nc) + w * term(ne)**2
              ENDDO
            ENDDO
          ENDDO

          nc = lmax*ntype+1
          w = 0.5*wparint(neig,nk)
          DO ne = 1,ned
            term(ne) = energy(ne) - ener
            IF ( energy(ne).GT.ecmax(neig) ) term(ne) = ecmax(neig)-ener
            IF ( term(ne).lt.0.0e0 ) term(ne)=0.0e0
            g(ne,nc) = g(ne,nc) + w * term(ne)**2
          ENDDO

        ENDDO
      ENDDO

      END SUBROUTINE tetra_dos
      END MODULE m_tetrados
