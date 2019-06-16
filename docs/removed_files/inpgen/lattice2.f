      MODULE m_lattice
      use m_juDFT
!---------------------------------------------------------------------!
! eventually easy input of all 14 Bravais lattices in agreement with  !
! 'International Tables of Crystallography'                           !
! table 2.1.1, ( page13 in 3rd Ed. )                                  !
!---------------------------------------------------------------------!
      CONTAINS
      SUBROUTINE lattice2( 
     >                    buffer,xl_buffer,errfh,bfh,nline,
     <                    a1,a2,a3,aa,scale,mat,i_c,ios )

      USE m_constants
      IMPLICIT NONE

!==> Arguments
      INTEGER, INTENT (IN) :: errfh,bfh,nline
      INTEGER, INTENT (IN)                 :: xl_buffer
      CHARACTER(len=xl_buffer), INTENT(IN) :: buffer
      REAL,    INTENT (OUT) :: a1(3),a2(3),a3(3)
      REAL,    INTENT (OUT) :: aa                ! overall scaling constant
      REAL,    INTENT (OUT) :: scale(3),mat(3,3) ! for trigonal lattices
      INTEGER, INTENT (OUT) :: i_c,ios

!==> Local Variables
      CHARACTER(len=40) :: latsys
      REAL    :: a0,a_rho
      REAL    :: a,b,c
      REAL    :: alpha,beta,gamma
      REAL    :: c1(3),c2(3),c3(3)
      REAL    :: ar,br,cr,b1,b2,am(3,3)
      REAL    :: ca,cb,at
      INTEGER :: i,j,err,i1,i2
      LOGICAL :: noangles

      REAL, PARAMETER :: eps = 1.0e-7,
     &              sqrt2  =  1.4142135623730950,
     &              sqrt3  =  1.7320508075688773,
     &              sqrt32 =  0.86602540378444,
     &             msqrt32 = -0.86602540378444

      REAL :: lmat(3,3,8)
      DATA  lmat /  1.0,  0.0,  0.0,      ! 1: primitive     : P
     &              0.0,  1.0,  0.0,
     &              0.0,  0.0,  1.0,  
     +              0.0,  0.5,  0.5,      ! 2: face-centered : F
     &              0.5,  0.0,  0.5,
     &              0.5,  0.5,  0.0,
     +             -0.5,  0.5,  0.5,      ! 3: body-centered : I
     &              0.5, -0.5,  0.5,
     &              0.5,  0.5, -0.5,
     +            0.5, msqrt32, 0.0,      ! 4: hexagonal-P   : hP, hcp
     &            0.5,  sqrt32, 0.0,
     &            0.0,     0.0, 1.0,   
     +              0.0, -1.0,  1.0,      ! 5: hexagonal-R   : hR, trigonal
     &           sqrt32,  0.5,  1.0,
     &          msqrt32,  0.5,  1.0, 
     +              0.5, -0.5,  0.0,      ! 6: base-centered: S (C)
     &              0.5,  0.5,  0.0,
     &              0.0,  0.0,  1.0,
     +              0.5,  0.0, -0.5,      ! 7: base-centered: B
     &              0.0,  1.0,  0.0,
     &              0.5,  0.0,  0.5,
     +              1.0,  0.0,  0.0,      ! 8: base-centered: A
     &              0.0,  0.5,  0.5,
     &              0.0, -0.5,  0.5/

!===> namelists
      NAMELIST /lattice/ latsys,a0,a,b,c,alpha,beta,gamma

      noangles = .false.
      latsys = ' ' ; a0 = 0.0   
      a = 0.0      ; b = 0.0    ; c = 0.0 
      alpha = 0.0  ; beta = 0.0 ; gamma = 0.0   
      scale = 0.0  ; mat = 0.0
 
      READ (bfh,lattice,err=911,end=911,iostat=ios)
      
      IF ( abs(a0) < eps ) a0 = 1.0
      IF ( abs(a)  < eps ) a  = 1.0
      IF ( abs(b)  < eps ) b  = a
      IF ( abs(c)  < eps ) c  = a
      IF ( abs(alpha)  < eps ) alpha  = 90.0
      IF ( abs(beta)   < eps ) beta   = 90.0
      IF ( abs(gamma)  < eps ) gamma  = 90.0

      IF ( alpha > pi_const ) THEN     ! deg
        ar = alpha * pi_const / 180.0 
        br = beta  * pi_const / 180.0 
        cr = gamma * pi_const / 180.0 
      ELSE                       ! radians
        ar = alpha
        br = beta
        cr = gamma
      ENDIF

      scale(1) = a
      scale(2) = b
      scale(3) = c

      latsys = ADJUSTL(latsys)

!===>  1: cubic-P          (cP) sc

      IF ( latsys =='cubic-P'.OR.latsys =='cP'.OR.latsys =='sc'.OR.
     &     latsys =='simple-cubic' ) THEN

        noangles=.true.
        i_c = 1
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0

        IF ( a.NE.b .OR. a.NE.c ) err = 11
        IF ( ar.NE.br .OR. ar.NE.cr .OR. ar.NE.(pi_const/2.0) ) err = 12

!===>  2: cubic-F          (cF) fcc

      ELSEIF ( latsys =='cubic-F'.OR.latsys =='cF'.OR.latsys =='fcc'.OR.
     &         latsys =='face-centered-cubic' ) THEN

        noangles=.true.
        i_c = 2
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0

        IF ( a.NE.b .OR. a.NE.c ) err = 21

!===>  3: cubic-I          (cI) bcc

      ELSEIF ( latsys =='cubic-I'.OR.latsys =='cI'.OR.latsys =='bcc'.OR.
     &         latsys =='body-centered-cubic' ) THEN

        noangles=.true.
        i_c = 3
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0

        IF ( a.NE.b .OR. a.NE.c ) err = 31

!===>  4: hexagonal-P      (hP) hcp

      ELSEIF ( latsys =='hexagonal-P'.OR.latsys =='hP'.OR.latsys =='hcp'
     &                               .OR.latsys =='hexagonal' ) THEN

        noangles=.true.
        i_c = 4
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0

        IF ( a.NE.b ) err = 41

!===>  4.1 : hexagonal-P   60 degrees variant

      ELSEIF ( latsys =='hdp' ) THEN

        noangles=.true.
        i_c = 4
        a1 =  lmat((/2,1,3/),1,i_c) * scale(:)
        a2 = -lmat((/2,1,3/),2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0
        i_c = 9

        IF ( a.NE.b ) err = 41

!===>  5.1: hexagonal-R      (hR) trigonal,( rhombohedral )

      ELSEIF ( latsys =='hexagonal-R'.OR.latsys =='hR'.OR.latsys =='r'.
     &      OR.latsys =='R'.OR.latsys =='rhombohedral'.OR.
     &         latsys =='rho'.OR.latsys =='trigonal' ) THEN

        noangles=.false.
        i_c = 5
        a1 = lmat(:,1,i_c)
        a2 = lmat(:,2,i_c)
        a3 = lmat(:,3,i_c)
        
        IF ( a.NE.b ) err = 51
        IF ( b.NE.c ) err = 51
        IF ( alpha.EQ.0.0  .OR. 
     &       alpha.NE.beta .OR. alpha.NE.gamma ) err = 52

        at = sqrt( 2.0 / 3.0 * ( 1 - cos(ar) ) )
        am(:,1) = a1 ; am(:,2) = a2 ; am(:,3) = a3
        am(1,:) = at * am(1,:)
        am(2,:) = at * am(2,:)
        am(3,:) = cos(asin(at)) * am(3,:)
        a1 = am(:,1)*a ; a2 = am(:,2)*b ; a3 = am(:,3)*c
        scale = 1.0

        CALL angles( am )
 
!===>  5.2: hexagonal-R      (hR) trigonal,( rhombohedral )

      ELSEIF (latsys =='hexagonal-R2'.OR.latsys =='hR2'.OR.latsys =='r2'
     &    .OR.latsys =='R2'.OR.latsys =='rhombohedral2'.OR.
     &        latsys =='trigonal2' ) THEN

        noangles=.false.
        i_c = 5
        a1 = lmat(:,1,i_c)
        a2 = lmat(:,2,i_c)
        a3 = lmat(:,3,i_c)

        IF ( a.NE.b ) err = 53
        IF ( alpha.NE.90 ) err = 54
        IF ( alpha.NE.beta .OR. gamma.NE.120 ) err = 54

        mat(:,1) = a*lmat(:,1,4) ! to transfer atom
        mat(:,2) = a*lmat(:,2,4) ! positions hex -> trig
        mat(:,3) = c*lmat(:,3,4) ! in struct_inp.f

! transform hex -> rho
        a_rho = sqrt( 3*a*a + c*c)/3.0
        ar = acos( 1. - 9./(6.+2*(c/a)**2) )
        at = sqrt( 2.0 / 3.0 * ( 1 - cos(ar) ) )

        a = a_rho * at
        b = a_rho * at
        c = a_rho * cos( asin(at) )

        am(:,1) = a1 ; am(:,2) = a2 ; am(:,3) = a3
        am(1,:) = a*am(1,:)
        am(2,:) = b*am(2,:)
        am(3,:) = c*am(3,:)
        a1 = am(:,1) ; a2 = am(:,2) ; a3 = am(:,3)
        scale = 1.0

        CALL angles( am )

!===>  6: tetragonal-P     (tP) st

      ELSEIF ( latsys =='tetragonal-P'.OR.latsys =='st'.OR.latsys =='tP'
     &     .OR.latsys =='simple-tetragonal' ) THEN

        noangles=.true.
        i_c = 1
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0
          
        IF ( a.NE.b ) err = 61
        IF ( ar.NE.br .OR. ar.NE.cr .OR. ar.NE.(pi_const/2.0)  ) err= 62

!===>  7: tetragonal-I     (tI) bct

      ELSEIF (latsys =='tetragonal-I'.OR.latsys =='tI'.OR.latsys =='bct'
     &    .OR.latsys =='body-centered-tetragonal' ) THEN

        noangles=.true.
        i_c = 3
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0
          
        IF ( a.NE.b ) err = 61

!===>  8: orthorhombic-P   (oP) 

      ELSEIF ( latsys =='orthorhombic-P'.OR.latsys =='oP'.OR.
     &         latsys =='simple-orthorhombic' ) THEN

        noangles=.true.
        i_c = 1
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0
          
!===>  9: orthorhombic-F   (oF) 

      ELSEIF ( latsys =='orthorhombic-F'.OR.latsys =='oF'.OR.
     &         latsys =='orF'.OR.
     &         latsys =='face-centered-orthorhombic' ) THEN

        noangles=.true.
        i_c = 2
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0
          
!===> 10: orthorhombic-I   (oI) 

      ELSEIF ( latsys =='orthorhombic-I'.OR.latsys =='oI'.OR.
     &         latsys =='orI'.OR.
     &         latsys =='body-centered-orthorhombic' ) THEN

        noangles=.true.
        i_c = 3
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0
          
!===> 11: orthorhombic-S   (oS) (oC) 

      ELSEIF ( latsys =='orthorhombic-S'.OR.latsys =='orthorhombic-C'.or
     &        .latsys =='oS'.OR.latsys =='oC'.OR.latsys =='orC'.OR.
     &         latsys =='base-centered-orthorhombic' ) THEN

        noangles=.true.
        i_c = 6
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0
          
!===> 11a: orthorhombic-A   (oA)

      ELSEIF ( latsys =='orthorhombic-A'.OR.latsys =='oA'.OR
     &        .latsys =='orA'.OR.
     &         latsys =='base-centered-orthorhombic2' ) THEN

        noangles=.true.
        i_c = 8
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0

!===> 11b: orthorhombic-B   (oB)

      ELSEIF ( latsys =='orthorhombic-B'.OR.latsys =='oB'.OR
     &        .latsys =='orB'.OR.
     &         latsys =='base-centered-orthorhombic3' ) THEN

        noangles=.true.
        i_c = 7
        a1 = lmat(:,1,i_c) * scale(:)
        a2 = lmat(:,2,i_c) * scale(:)
        a3 = lmat(:,3,i_c) * scale(:)
        scale = 1.0

!===> 12: monoclinic-P     (mP) 
      ELSEIF ( latsys =='monoclinic-P'.OR.latsys =='mP'.OR
     &        .latsys =='moP'.OR.
     &         latsys =='simple-monoclinic' ) THEN

        noangles=.false.
        i_c = 1

        IF ( (abs(alpha-90.0)<eps).AND.(abs(beta-90.0)<eps) ) THEN
           IF ( ABS(gamma - 90.0) <eps )  CALL juDFT_error
     +          ("no monoclinic angle!",calledby ="lattice2")
        ELSE 
           CALL juDFT_error
     +          ("lattice2: Please take gamma as monoclinic angle!"
     +          ,calledby ="lattice2")
        ENDIF  
        CALL brvmat ( alpha, beta, gamma, am )
        a1 = a*am(:,1)
        a2 = b*am(:,2)
        a3 = c*am(:,3) 
        scale = 1.0

        CALL angles( am )

!===> 13: monoclinic-C (mC) 
      ELSEIF ( latsys =='monoclinic-C'.OR.latsys =='mC'.OR
     &        .latsys =='moC'.OR.
     &         latsys =='centered-monoclinic' ) THEN

         CALL juDFT_error
     +    ("Monoclinic setup",hint
     +        ="Please use gamma as monoclinic angle"//
     +        " and center on A or B !",calledby ="lattice2")

!===> 13a monoclinic-A (mA)
      ELSEIF ( latsys =='monoclinic-A'.OR.latsys =='mA'.OR
     &        .latsys =='moA'.OR.
     &         latsys =='centered-monoclinic2' ) THEN

        noangles=.false.
        IF ( (abs(alpha-90.0)<eps).AND.(abs(beta-90.0)<eps) ) THEN
           IF ( abs(gamma - 90.0) <eps )  CALL juDFT_error
     +          ("no monoclinic angle!",calledby ="lattice2")
        ELSE
           CALL juDFT_error("Please take gamma as monoclinic angle!"
     +          ,calledby ="lattice2")
        ENDIF
        CALL brvmat ( alpha, beta, gamma, am )
        am(:,1) = a * am(:,1)
        am(:,2) = b * am(:,2)
        am(:,3) = c * am(:,3)
        i_c = 8
        am = matmul ( am, lmat(:,:,i_c) )
        a1 = am(:,1)
        a2 = am(:,2)
        a3 = am(:,3)
        CALL angles( am )
        scale = 1.0 

!===> 13b monoclinic-B (mB)
      ELSEIF ( latsys =='monoclinic-B'.OR.latsys =='mB'.OR
     &        .latsys =='moB'.OR.
     &         latsys =='centered-monoclinic3' ) THEN

        noangles=.false.
        IF ( (abs(alpha-90.0)<eps).AND.(abs(beta-90.0)<eps) ) THEN
           IF ( ABS(gamma - 90.0) <eps )  CALL juDFT_error
     +          ("no monoclinic angle!",calledby ="lattice2")
        ELSE
           CALL juDFT_error("Please take gamma as monoclinic angle!"
     +          ,calledby ="lattice2")
        ENDIF
        CALL brvmat ( alpha, beta, gamma, am )
        am(:,1) = a * am(:,1)
        am(:,2) = b * am(:,2)
        am(:,3) = c * am(:,3)
        i_c = 7
        am = matmul ( am, lmat(:,:,i_c) )
        a1 = am(:,1)
        a2 = am(:,2)
        a3 = am(:,3)
        CALL angles( am )
        scale = 1.0 

!===> 13c monoclinic-I (mI)
      ELSEIF ( latsys =='monoclinic-I'.OR.latsys =='mI'.OR
     &        .latsys =='moI' ) THEN

        noangles=.false.
        IF ( (abs(alpha-90.0)<eps).AND.(abs(beta-90.0)<eps) ) THEN
           IF ( ABS(gamma - 90.0) <eps )  CALL juDFT_error
     +          ("no monoclinic angle!",calledby ="lattice2")
        ELSE
           CALL juDFT_error("Please take gamma as monoclinic angle!"
     +          ,calledby ="lattice2")
        ENDIF
        CALL brvmat ( alpha, beta, gamma, am )
        am(:,1) = a * am(:,1)
        am(:,2) = b * am(:,2)
        am(:,3) = c * am(:,3)
        i_c = 3
        am = matmul ( am, lmat(:,:,i_c) )
        a1 = am(:,1)
        a2 = am(:,2)
        a3 = am(:,3)
        CALL angles( am )
        scale = 1.0 

!===> 14: triclinic        (aP) 

      ELSEIF ( latsys =='aP' .OR. latsys =='triclinic' .OR.
     &         latsys =='tcl' )  THEN

        noangles=.false.
        i_c = 1
        CALL brvmat ( alpha, beta, gamma, am )
        a1 = a*am(:,1)
        a2 = b*am(:,2)
        a3 = c*am(:,3)
        CALL angles( am )
        scale = 1.0
          
      ELSE
          WRITE (errfh,*)
          WRITE (errfh,*) '*** unknown lattice system latsys=',latsys,
     &                    'in line',nline,'.'
          WRITE (errfh,*) 
     &    '***  No reason for panic, it is probably because'
          WRITE (errfh,*) 
     & '***  subroutine lattice2 in struct_input.f is unfinished. ***'
        ios = 1
        RETURN
      ENDIF

      IF ( noangles .AND.
     &    ( abs(alpha - 90.0 ) > eps .OR.
     &      abs(beta  - 90.0 ) > eps .OR.
     &      abs(gamma - 90.0 ) > eps )    ) THEN
        WRITE (errfh,*)
        WRITE (errfh,*) 'ERROR in &lattice ... /. ',
     & 'For the given lattice system all angles should be 90deg.'
        WRITE (errfh,*)
      ENDIF

      IF (abs(scale(1)) < eps) THEN
        CALL juDFT_error("Illegal program section reached!"
     +          ,calledby ="lattice2")
        scale(1) = a
        scale(2) = b
        scale(3) = c
      ENDIF
      aa = a0

 911  CONTINUE

      RETURN
      END SUBROUTINE lattice2
!
!**********************************************************************!
!
      SUBROUTINE angles( am )
!----------------------------------------------------------------------!
!     given an bravais matrix (am), calculate length of basis vectors  !
!     and angles beween them ; write results to standard output        !
!----------------------------------------------------------------------!

      USE m_constants
      IMPLICIT NONE

      REAL, INTENT(IN) :: am(3,3)

      REAL     ca,al(3)
      INTEGER  i,j,i1,i2

      al = 0.0
      DO j = 1, 3
        DO i = 1, 3
          al(j) = al(j) + am(i,j)*am(i,j)
        ENDDO
        al(j) = sqrt(al(j))
      ENDDO

      DO j = 1, 3
        WRITE (6,'("vector ",i1," : ",3f9.5,5x," length : ",f9.5)') 
     &                                              j,am(:,j),al(j)
      ENDDO

      DO i1 = 1, 2 
        DO i2 = i1+1, 3
          ca = 0.0
          DO i = 1, 3
            ca = ca + am(i,i1)*am(i,i2)
          ENDDO
          ca  = ca/(al(i1)*al(i2))
          ca = acos(ca)*180/pi_const
    
          WRITE (6,'("angle between vectors (",i1,",",i1,") =",f9.5)') 
     &                                                        i1,i2,ca
        ENDDO
      ENDDO

      END SUBROUTINE angles
!
!**********************************************************************!
!
      SUBROUTINE brvmat ( alpha, beta, gamma, am )
!----------------------------------------------------------------------!
!     given the angles alpha, beta and gamma, set up an matrix 'am'    !
!     with 3 unit vectors and these angles between them. The first     !
!     unit vector poins in (1,0,0) direction                gb`05      !
!----------------------------------------------------------------------!

      USE m_constants
      IMPLICIT NONE

      REAL, INTENT (IN)  :: alpha, beta, gamma
      REAL, INTENT (OUT) :: am(3,3)

      REAL ca,cb,cg,sg,c1,c2

      ca = cos(alpha*pi_const/180); cb = cos(beta*pi_const/180)
      cg = cos(gamma*pi_const/180); sg = sin(gamma*pi_const/180)
      c1 = (ca - cg*cb ) / sg
      c2 = sqrt( 1 - cb**2 - c1**2 ) 

      am(1,1) = 1.0 ; am(2,1) = 0.0 ; am(3,1) = 0.0
      am(1,2) = cg  ; am(2,2) = sg  ; am(3,2) = 0.0
      am(1,3) = cb  ; am(2,3) = c1  ; am(3,3) = c2

      END SUBROUTINE brvmat

      END MODULE m_lattice
