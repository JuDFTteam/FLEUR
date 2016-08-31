MODULE m_anglso
contains
  COMPLEX FUNCTION anglso(theta,phi,l1,m1,is1,l2,m2,is2)
    USE m_juDFT
    !
    ! calculates spin-orbit matrix for theta,phi =/= 0
    !
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: is1,is2,l1,l2,m1,m2
    REAL,    INTENT(IN) :: theta,phi             
    !     ..
    !     .. Local Scalars ..
    REAL sgm1,sgm2,xlz,xlpl,xlmn,angl_r,angl_i
    LOGICAL :: l_standard_euler_angles
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC abs,REAL,sqrt,isign
    !
    anglso = CMPLX(0.0,0.0)
    IF (l1.NE.l2) THEN
       RETURN
    ENDIF
    !

    l_standard_euler_angles=.FALSE.	

    sgm1 = is1
    sgm2 = is2
    IF (l1.LT.0) THEN
       WRITE (6,FMT=*)&
            &      ' PROGRAM STOPS IN ANGLSO ( L < 0 ) .'
       WRITE (6,FMT=*) ' L1 =',l1,'    L2 =',l2
       CALL juDFT_error("ANGLSO (L <0 )",calledby="anglso")
    ELSE IF ((ABS(m1).GT.l1) .OR. (ABS(m2).GT.l2)) THEN
       WRITE (6,FMT=*)&
            &      ' PROGRAM STOPS IN ANGLSO ( M < L OR L < M )'
       WRITE (6,FMT=*) ' L1 =',l1,'    L2 =',l2
       WRITE (6,FMT=*) ' M1 =',m1,'    M2 =',m2
       CALL juDFT_error("ANGLSO ( M < L OR L < M )",calledby="anglso"&
            &         )
    ELSE IF ((is1.NE.-1.AND.is1.NE.1) .OR.&
         &         (is2.NE.-1.AND.is2.NE.1)) THEN
       WRITE (6,FMT=*)&
            &      ' PROGRAM STOPS IN ANGLSO ( S >< +-1/2 ) .'
       WRITE (6,FMT=*) ' S1 =',0.5*sgm1,'    S2 =',0.5*sgm2
       CALL juDFT_error("ANGLSO ( S >< +-1/2 )",calledby ="anglso")
    END IF
    !
    ! lz,l+,l-
    !
    xlz = 0.0
    xlpl= 0.0
    xlmn= 0.0
    ! is1.eq.is2-2 -> <-| |+> => l+   
    IF (m1.EQ.m2+1) THEN
       xlpl = SQRT(REAL((l2-m2)* (l2+m2+1)))
       ! is1.eq.is2+2 -> <+| |-> => l-   
    ELSE IF (m1.EQ.m2-1) THEN
       xlmn = SQRT(REAL((l2+m2)* (l2-m2+1)))
       ! is1.eq.is2 -> <+| |+> => lz   
    ELSE IF (m1.EQ.m2  ) THEN
       xlz  = m2
    END IF
    !
    ! rotated spin-orbit angular matrix
    ! <1| |1> or <2| |2>          
    !

    !     If l_standard_euler_angles is set to .false., the old version
    !     of this subroutine is used, which works fine if only MAE is
    !     needed. In the old version, theta and phi specify the direction
    !     of spin quantization. However, there remains the degree of freedom
    !     of rotation around the spin quantization direction. The new version
    !     of this subroutine, which is switched on by setting l_standard_euler_angles
    !     to true, rotates x, y and z axes first by phi around z counterclockwisely
    !     and then the resulting axes x', y' and z' by theta around y' 
    !     counterclockwisely to obtain the final coordinate frame x'', y'' and z''.
    !     This new version is useful when one is interested also in the spin
    !     perpendicular to the spin quantization axis. The old version seems to be
    !     different from the new one as follows when theta and phi differ from zero:
    !     x''_old = -x''_new
    !     y''_old = -y''_new
    !     z''_old = z''
    !     If theta and phi are zero, there is no difference, leading thus to a
    !     discontinuous Euler transformation in the old version and a continuous
    !     transformation in the new version.
    !     Both versions do not differ for MAE, but if the standard Pauli matrix
    !     is used to obtain x or y components of spin, the spin is discontinuous
    !     in the old version as a function of theta and phi.

    IF(l_standard_euler_angles)THEN	

       IF (is1.EQ.is2) THEN
          angl_r = isign(1,is1) * ( COS(theta)*xlz +&
               &                          0.5*SIN(theta)*COS(phi)*(xlmn + xlpl) )
          angl_i = isign(1,is1)*0.5*SIN(theta)*SIN(phi)*(xlmn - xlpl)
          ! <1| |2>
       ELSEIF (is1.EQ.is2+2) THEN
          angl_r =  - SIN(theta)*xlz +  COS(phi)*(&
               &             COS(theta/2.)**2 * xlmn - SIN(theta/2.)**2 * xlpl )
          angl_i = SIN(phi)*( &
               &              COS(theta/2.)**2 * xlmn + SIN(theta/2.)**2 * xlpl )
          ! <2| |1>
       ELSEIF (is1.EQ.is2-2) THEN
          angl_r = - SIN(theta)*xlz +  COS(phi)*(&
               &            + COS(theta/2.)**2 * xlpl - SIN(theta/2.)**2 * xlmn )
          angl_i =  SIN(phi)*( &
               &             - COS(theta/2.)**2 * xlpl - SIN(theta/2.)**2 * xlmn )
       ELSE
          angl_r = 0.0
          angl_i = 0.0
       ENDIF

    ELSE 	
       IF (is1.EQ.is2) THEN
          angl_r = isign(1,is1) * ( COS(theta)*xlz +&
               &                          0.5*SIN(theta)*COS(phi)*(xlmn + xlpl) )
          angl_i = isign(1,is1)*0.5*SIN(theta)*SIN(phi)*(xlmn - xlpl)
          ! <1| |2>
       ELSEIF (is1.EQ.is2+2) THEN
          angl_r =  SIN(theta)*xlz +  COS(phi)*(&
               &            - COS(theta/2.)**2 * xlmn + SIN(theta/2.)**2 * xlpl )
          angl_i = -SIN(phi)*( &
               &              COS(theta/2.)**2 * xlmn + SIN(theta/2.)**2 * xlpl )
          ! <2| |1>
       ELSEIF (is1.EQ.is2-2) THEN
          angl_r =  SIN(theta)*xlz +  COS(phi)*(&
               &            - COS(theta/2.)**2 * xlpl + SIN(theta/2.)**2 * xlmn )
          angl_i =  SIN(phi)*( &
               &              COS(theta/2.)**2 * xlpl + SIN(theta/2.)**2 * xlmn )
       ELSE
          angl_r = 0.0
          angl_i = 0.0
       ENDIF
    ENDIF
    !
    anglso = CMPLX(angl_r,angl_i)

    RETURN
  END FUNCTION anglso
END MODULE m_anglso
