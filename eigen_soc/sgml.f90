MODULE m_sgml
  CONTAINS
  REAL FUNCTION sgml(l1,m1,is1,l2,m2,is2)
    USE m_juDFT
    !
    ! FUNCTION SGML ******************************************************
    !
    !      THIS FUNCTION CALCULATES ANGULAR PART OF THE MATRIX ELEMENT
    !
    !      ( < L1 MP1 IS1 ! SGM*L ! L2 MP2 IS2 > ) OF SPIN-ORBIT COUPLING
    !
    !      FOR COMPLEX :: SPHERICAL HARMONICS FOR :
    !      L >= 0 ,-L <= M <= L , SGM = 1 OR -1 , IS = 1   OR -1
    !
    !                             SATOSHI TAKIZAWA, ISSP, MAR 1990
    !            MODIFIED BY      STEFAN BL"UGEL  , ISSP, MAR 1990
    ! *********************************************************************
    !
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    INTEGER is1,is2,l1,l2,m1,m2
    !     ..
    !     .. Local Scalars ..
    REAL sgm1,sgm2
    !     ..
    !     ..
    IF (l1.NE.l2) THEN
       sgml = 0.0
       RETURN
    ELSE
       sgm1 = is1
       sgm2 = is2
       IF (l1.LT.0) THEN
          WRITE (6,FMT=*) ' PROGRAM STOPS IN FUNCTION SGML ( L < 0 ) .'
          WRITE (6,FMT=*) ' L1 =',l1,'    L2 =',l2
          CALL juDFT_error("SGMLR",calledby="sgml")
       ELSE IF ((ABS(m1).GT.l1) .OR. (ABS(m2).GT.l2)) THEN
          WRITE (6,FMT=*) ' PROGRAM STOPS IN SGMLC ( jij%M < L OR L < jij%M )'
          WRITE (6,FMT=*) ' L1 =',l1,'    L2 =',l2
          WRITE (6,FMT=*) ' M1 =',m1,'    M2 =',m2
          CALL juDFT_error("SGML",calledby="sgml")
       ELSE IF ((is1.NE.-1.AND.is1.NE.1) .OR. (is2.NE.-1.AND.is2.NE.1)) THEN
          WRITE (6,FMT=*) ' PROGRAM STOPS IN FUNCTION SGMLC ( S >< +-1/2 ) .'
          WRITE (6,FMT=*) ' S1 =',0.5*sgm1,'    S2 =',0.5*sgm2
          CALL juDFT_error("SGML",calledby="sgml")
       END IF
       !
       !
       IF (m1.EQ.m2+1 .AND. is1.EQ.is2-2) THEN
          sgml = SQRT(REAL((l2-m2)* (l2+m2+1)))
       ELSE IF (m1.EQ.m2-1 .AND. is1.EQ.is2+2) THEN
          sgml = SQRT(REAL((l2+m2)* (l2-m2+1)))
       ELSE IF (m1.EQ.m2 .AND. is1.EQ.is2) THEN
          sgml = m2*sgm2
       ELSE
          sgml = 0.0
       END IF
       RETURN
    END IF
    !
    RETURN
  END FUNCTION sgml
END MODULE m_sgml
