!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_bravais
      use m_juDFT
!----------------------------------------------------------------------!
! given a Bravais-matrix amat, determine the lattice system and type   !
! (idsyst,idtype)                                              gb`05   !
!----------------------------------------------------------------------!
      CONTAINS
      SUBROUTINE bravais(
     >                   amat,
     <                   idsyst,idtype)

      IMPLICIT NONE

      REAL,    INTENT (IN)  :: amat(3,3)
      INTEGER, INTENT (OUT) :: idsyst,idtype

      REAL a(3),b(3),c(3),g(3),e(3),f(3)
      REAL sa, sb, sc, al, be, ga
      REAL, PARAMETER :: eps = 1.0e-5
      REAL, PARAMETER :: thrd = -1./3.
      INTEGER i,j,k
      LOGICAL l_ab, l_bc, l_ac, al_be, be_ga, al_ga

      CHARACTER (len=12) :: c_sy(7)
      CHARACTER (len=15) :: c_ty(6)
      c_sy = (/'cubic       ','tetragonal  ','orthorhombic',
     +         'hexagonal   ','trigonal    ','monoclinic  ',
     +         'triclinic   '/)
      c_ty = (/'primitive      ','body centered  ','face centered  ',
     +         'A-face centered','B-face centered','C-face centered'/)

      a(:) = amat(:,1) ; b(:) = amat(:,2) ; c(:) = amat(:,3)
      sa = SQRT( DOT_PRODUCT( a, a) )
      sb = SQRT( DOT_PRODUCT( b, b) )
      sc = SQRT( DOT_PRODUCT( c, c) ) 
      al = DOT_PRODUCT( b, c) / ( sb * sc )
      be = DOT_PRODUCT( a, c) / ( sa * sc )
      ga = DOT_PRODUCT( b, a) / ( sb * sa )
      write (*,*) sa,sb,sc,al,be,ga

      l_ab = .false. ; l_bc = .false. ; l_ac = .false.
      al_be = .false. ; be_ga = .false. ; al_ga = .false.
      IF ( ABS(sa-sb) < eps ) l_ab = .true.
      IF ( ABS(sb-sc) < eps ) l_bc = .true.
      IF ( ABS(sa-sc) < eps ) l_ac = .true.
      IF ( ABS(al-be) < eps ) al_be = .true.
      IF ( ABS(be-ga) < eps ) be_ga = .true.
      IF ( ABS(al-ga) < eps ) al_ga = .true.
      
      idsyst = 99 ; idtype = 99

      IF (l_ab.AND.l_bc)  THEN                ! all sides equal
        IF (al_be.AND.be_ga) THEN             ! all angles equal
          IF ( ABS(al) < eps )  THEN          ! alpha = 90 deg
            idsyst = 1 ; idtype = 1           !       --> simple cubic
          ELSEIF ( ABS(al-thrd) < eps ) THEN  ! alpha = 109 deg
            idsyst = 1 ; idtype = 2           !       --> bcc
          ELSEIF ( ABS(al-0.5) < eps ) THEN   ! alpha = 60  deg
            idsyst = 1 ; idtype = 3           !       --> fcc
          ELSE
            idsyst = 5 ; idtype = 1           ! -->  trigonal (rhomboedric) 
          ENDIF
        ELSEIF (al_be.OR.be_ga.OR.al_ga) THEN ! two angles equal
          idsyst = 2 ; idtype = 2             !  -->  tetragonal - I
        ELSE
          idsyst = 3 ; idtype = 2             !  -->  orthorhombic - I
        ENDIF
      ENDIF

      IF (idsyst == 99) THEN                  ! continue the search

      IF (l_ab.OR.l_bc.OR.l_ac)  THEN         ! two sides equal
                                              ! hexagonal or tetragonal or base-centered
        IF (al_be.AND.be_ga) THEN             ! all angles equal
          IF ( ABS(al) < eps )  THEN          ! alpha = 90  deg
            idsyst = 2 ; idtype = 1           !       --> simple tetragonal
          ELSE
            IF (l_ab) THEN
             idsyst = 6 ; idtype = 6          ! special monoclinic - C
            ELSEIF (l_bc) THEN
             idsyst = 6 ; idtype = 4          ! special monoclinic - A
            ELSEIF (l_ac) THEN
             idsyst = 6 ; idtype = 5          ! special monoclinic - B
            ENDIF
          ENDIF
        ELSEIF (.NOT.(al_be.OR.be_ga.OR.al_ga)) THEN  ! all angles different
          idsyst = 7 ; idtype = 1                     ! triclinic
        ELSE
          IF ( al_be.AND.( ABS(ABS(ga)-0.5) < eps ) ) THEN
            idsyst = 4 ; idtype = 1 ! hexagonal
          ELSEIF ((ABS(al)<eps).OR.(ABS(be)<eps).OR.(ABS(ga)<eps)) THEN  ! one is 90 deg
            IF ( (ABS(al)<eps).AND.(ABS(be)<eps) ) THEN
              idsyst = 3 ; idtype = 6                         ! orthorhombic - C
            ELSEIF ( (ABS(ga)<eps).AND.(ABS(be)<eps) ) THEN
              idsyst = 3 ; idtype = 4                         ! orthorhombic - A
            ELSEIF ( (ABS(al)<eps).AND.(ABS(ga)<eps) ) THEN
              idsyst = 3 ; idtype = 5                         ! orthorhombic - B
            ELSE
              idsyst = 6 ; idtype = 1            ! simple monoclinic
            ENDIF
          ELSE                                                           ! none is 90 deg
            IF (al_be) THEN
             idsyst = 6 ; idtype = 6             ! monoclinic - C
            ELSEIF (be_ga) THEN
             idsyst = 6 ; idtype = 4             ! monoclinic - A
            ELSEIF (al_ga) THEN
             idsyst = 6 ; idtype = 5             ! monoclinic - B
            ELSE
              idsyst = 6 ; idtype = 1            ! simple monoclinic
            ENDIF
          ENDIF
        ENDIF

      ELSE                                    ! orthorhombic or tricinic or monoclinic

        IF (al_be.AND.be_ga) THEN             ! all angles equal
          IF ( ABS(al) < eps )  THEN          ! angles 90 deg
           idsyst = 3 ; idtype = 1            ! --> simple orthorhombic
          ELSE
           idsyst = 7 ; idtype = 1            ! triclinic
          ENDIF
        ELSEIF (.NOT.(al_be.OR.be_ga.OR.al_ga)) THEN ! all angles different
          e = a + b - c ; f = b + c - a ; g = a + c - b
          IF ( (DOT_PRODUCT( e, f) == 0).AND.
     +         (DOT_PRODUCT( e, g) == 0).AND.
     +         (DOT_PRODUCT( g, f) == 0) )  THEN
             idsyst = 3 ; idtype = 3          ! --> face-centered orthorhombic
          ELSE
             idsyst = 7 ; idtype = 1         ! triclinic
          ENDIF
        ELSE
          idsyst = 6 ; idtype = 1   ! simple monoclinic 
        ENDIF

      ENDIF

      ENDIF
     
      IF ((idsyst == 99).OR.(idtype == 99) ) CALL juDFT_error("bravais!"
     +     ,calledby ="bravais")
 10   WRITE(*,*) c_ty(idtype),' ',c_sy(idsyst)

      END SUBROUTINE bravais
      END MODULE m_bravais
