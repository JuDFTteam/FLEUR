!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fertri
   !
   !     calculates fermi energy and weights using triangular method
   !
   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_triang
   USE m_maketetra
   USE m_tetraef
   USE m_dosef
   USE m_dosint
   USE m_doswt

   IMPLICIT NONE
   CONTAINS
   SUBROUTINE fertri(input,kpts,irank,ne,nkpt,jspins,zc,eig,bk,sfac,&
                     ef,seigv,w)

      TYPE(t_input),INTENT(IN):: input
      TYPE(t_kpts), INTENT(IN):: kpts
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)    :: nkpt,jspins,irank
      REAL,    INTENT (IN)    :: zc,sfac
      REAL,    INTENT (OUT)   :: seigv
      REAL,    INTENT (INOUT) :: ef
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN)    :: ne(:,:)!(nkptd,jspd)
      REAL,    INTENT (IN)    :: bk(:,:) !(3,nkptd)
      REAL,    INTENT (OUT)   :: w(:,:,:) !(neigd,nkptd,jspd)
      REAL,    INTENT (INOUT) :: eig(:,:,:)!(neigd,nkptd,jspd)
!     ..
!     .. Local Scalars ..
      REAL chmom,ct,de,del,dez,ei,emax,emin,s,s1,workf
      REAL lb,ub,e_set
      LOGICAL film
      INTEGER i,ic,j,jsp,k,neig
      INTEGER ntria,ntetra      ! number of triangles & tetrahedrons
      REAL as                   ! total area covered by triangles
!     ..
!     .. Local Arrays ..
      INTEGER itria(3,2*size(w,2))  ! index of k-points that are corner points of a triangle
      REAL    atr(2*size(w,2))      ! area of a triangle
      INTEGER itetra(4,6*size(w,2)) ! ditto for tetrahedrons
      REAL    voltet(6*nkpt)
      INTEGER nemax(2)
!     ..
!     .. Data statements ..
      DATA de/5.0e-3/

      IF ( irank == 0 ) THEN
        WRITE (oUnit,FMT=8000)
      END IF
 8000 FORMAT (/,/,10x,'linear triangular method')

      film = .true.
      CALL triang(bk,nkpt,itria,ntria,atr,as,film)!keep
!
!--->   clear w and set eig=-9999.9
      e_set = -9999.9
      IF (.NOT.film) e_set = 1.0e10
      DO jsp = 1,jspins
         nemax(jsp) = 0.0
         DO k = 1,nkpt
            nemax(jsp) = max0(nemax(jsp),ne(k,jsp))
            DO i = 1,ne(k,jsp)
               w(i,k,jsp) = 0.
            ENDDO
            DO i = ne(k,jsp)+1,size(w,1)
               w(i,k,jsp) = 0.
               eig(i,k,jsp) = e_set
            ENDDO
         ENDDO
      ENDDO
!
!      sfac = 2.0/real(jspins)
!
!--->   write results of triang

      IF (.not.film) THEN
        ntetra = kpts%ntet
        DO j = 1, ntetra
          itetra(1:4,j) = kpts%ntetra(1:4,j)
          voltet(j) = kpts%voltet(j) / ntetra
        END DO
        lb = MINVAL(eig(:,:,:)) - 0.01
        ub = ef + 0.2
        CALL tetra_ef(jspins,nkpt,lb,ub,eig,zc,sfac,ntetra,itetra,voltet,ef,w)
      ELSE

        DO i = 1,ntria
           atr(i) = atr(i)/as
        ENDDO
        IF ( irank == 0 ) THEN
          WRITE (oUnit,FMT=8010) ntria,as
          DO i = 1,ntria
            WRITE (oUnit,FMT=8020) i, (itria(j,i),j=1,3),atr(i)
          ENDDO
        END IF
 8010   FORMAT (/,10x,'triangular decomposition of brillouin zone:',/,&
                10x,'number of triangles=',i3,/,10x,&
                'total area of triangles=',f12.6,/,10x,&
                'no.,corners and (normalized) area of each triangle:',/)
 8020   FORMAT (10x,i3,3x,3i3,f14.6)
        IF ( irank == 0 ) THEN
          WRITE (oUnit,FMT=*) 'ef_hist=',ef
        END IF
        ei = ef
!jr     emin = -9999.9
        emin = +9999.9
        emax = -emin
        ic = 1
   90   IF (ic.GT.100) GO TO 230
        ic = ic + 1
!
!     results from triang are included here
!
        CALL dosint(ei,nemax,jspins,sfac,ntria,itria,atr,eig,ct)
!
        IF ( irank == 0 ) WRITE (oUnit,FMT=*) 'ct=',ct

        IF (ct.LT.zc) THEN            ! ei < ef
          emin = ei
          ei = ei + de
          IF (emin.GT.emax) GO TO 90
        ELSEIF (ct.GT.zc) THEN        ! ei > ef
          emax = ei
          ei = ei - de
          IF (emin.GT.emax) GO TO 90
        ENDIF
        IF (ct.NE.zc) THEN
          IF ( irank == 0 ) WRITE (oUnit,FMT=*) '2nd dosint'
!--->     refine ef to a value of 5 mry * (2**-20)
          iterate : DO i = 1, 40
            ei = 0.5* (emin+emax)
!
            CALL dosint(ei,nemax,jspins,sfac,ntria,itria,atr,eig,ct)
!
            IF ( irank == 0 ) WRITE (oUnit,FMT=*) 'i=',i,', ct=',ct
            IF ( ct == zc ) THEN
              EXIT iterate
            ELSEIF ( ct > zc ) THEN
              emax = ei
            ELSE
              emin = ei
            ENDIF
          ENDDO iterate
        ENDIF
        ef = ei
        del = emax - emin
        dez = zc - ct
        workf = -13.6058*2*ef
        IF ( irank == 0 ) THEN
          WRITE (oUnit,FMT=8030) ef,workf,del,dez
        END IF
 8030   FORMAT(/,10x,'fermi energy=',f10.5,' har',/,10x,'work function='&
               ,f10.5,' ev',/,10x,'uncertainity in energy and weights=',&
               2e16.6)
!
!--->   obtain dos at ef
!
        CALL dosef(ei,nemax,jspins,sfac,ntria,itria,atr,eig)
!
!--->   obtain weights needed for integration
!
        CALL doswt(ei,nemax,jspins,ntria,itria,atr,eig,w)

      ENDIF ! .NOT.film
!
!--->   write weights
!
!      DO 190 jsp = 1,jspins
!         neig = nemax(jsp)
!         DO 180 i = 1,neig
!            DO 170 k = 1,nkpt
!             WRITE (oUnit,FMT=*) 'w(',i,',',k,',',jsp,')=',w(i,k,jsp)
!  170       CONTINUE
!  180    CONTINUE
!  190 CONTINUE
!
!--->   obtain sum of weights and valence eigenvalues
!
      s1 = 0.
      seigv = 0.
      DO 220 jsp = 1,jspins
         s = 0.
         neig = nemax(jsp)
         DO 210 i = 1,neig
            DO 200 k = 1,nkpt
               s = s + w(i,k,jsp)
               seigv = seigv + w(i,k,jsp)*eig(i,k,jsp)
  200       CONTINUE
  210    CONTINUE
         s1 = s1 + s
  220 CONTINUE
      seigv = sfac*seigv
      chmom = s1 - jspins*s
      IF ( irank == 0 ) THEN
        WRITE (oUnit,FMT=8040) seigv,s1,chmom
      END IF
 8040 FORMAT (/,10x,'sum of valence eigenvalues=',f20.6,5x,&
             'sum of weights=',f10.6,/,10x,'moment=',f12.6)
      RETURN
!
  230 IF ( irank == 0 ) THEN
        WRITE (oUnit,FMT=8050) ei,ef,emin,emax,ct,zc
      END IF
 8050 FORMAT (/,/,10x,'error fertri: initial guess of ef off by 25 mry',&
             ' ei,ef,emin,emax,ct,zc',/,10x,6e16.7,/,10x,&
             'check number of bands')
      CALL juDFT_error("initial guess of ef off by 25 mry",calledby="fertri")
   END SUBROUTINE fertri
END MODULE m_fertri
