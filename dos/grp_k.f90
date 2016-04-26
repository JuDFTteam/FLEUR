MODULE m_grp_k
  USE m_juDFT

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: grp_k, euler

CONTAINS 

  SUBROUTINE grp_k(sym,mrot_k,cell,bk,nclass,nirr,char_table,&
       &    grpname,irrname,su,sp_alph,sp_beta)

    !************************************************************
    !
    !   Determines the group of k, returns the number of classes, the name of the
    !   group, the irreducible representations and the character table.
    !   All the groups are not implemented yet, and the identification is not tested for 
    !   the groups. 
    !
    !  The subroutine works also with double groups, however, no double groups are 
    !  tabulated at the moment.
    !
    ! 
    !   Character tables are taken from 
    !   T. Inui, Y Tanabe, and Y. Onodera, "Group theory and its applications in physics",
    !   Springer (1996)
    !
    !   Jussi Enkovaara, Juelich 2004
    !**************************************************************

    !      USE m_mrot2su
    USE m_inv3
    USE m_constants, ONLY : pi_const
    USE m_socsym,    ONLY : soc_sym, cross
    USE m_types
    IMPLICIT NONE

    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
   

    REAL,INTENT(IN)      :: bk(3)
    COMPLEX, INTENT(OUT) :: char_table(:,:)
    INTEGER, INTENT(OUT) :: mrot_k(:,:,:),nclass,nirr
    CHARACTER(LEN=7),  INTENT(OUT) :: grpname
    CHARACTER(LEN=5),  INTENT(OUT) :: irrname(:)
    COMPLEX, OPTIONAL, INTENT(OUT) :: su(:,:,:)
    REAL,    OPTIONAL, INTENT(IN)  :: sp_alph,sp_beta ! spin quant. axis

    ! locals
    INTEGER ::  nopk,c
    LOGICAL, ALLOCATABLE :: error(:)
    INTEGER :: mrot2(3,3,48)
    REAL    :: ktest(3),rtmp,rtmp2(4),kt(3)
    INTEGER :: mtest(3,3),mtmp(3,3),mtmpinv(3,3),munit(3,3)
    INTEGER :: n,n2 ,i,itmp
    INTEGER :: elem(48),belongs(48),members(48)
    LOGICAL :: soc,l_sorted
    ! COMPLEX :: sutmp(2,2), sutmpinv(2,2), su2(2,2,48)
    COMPLEX :: sutest(2,2)
    INTEGER :: d,det(48),rot(12),rot_type(48)
    REAL :: alpha,beta,gamma,theta(48),rax(4,48)
    REAL,PARAMETER :: eps=0.000001
    COMPLEX :: omega
    COMPLEX,PARAMETER :: one=CMPLX(1.0,0.0),zero=CMPLX(0.0,0.0)
    COMPLEX,PARAMETER :: imi=CMPLX(0.0,1.0)


    soc=.FALSE.
    IF (PRESENT(su)) soc=.TRUE.
    !sym%nop=SIZE(sym%mrot,3)

    ALLOCATE(error(sym%nop))
    error=.FALSE.
    ! Reduce the symmetry due to spin-orbit
    IF (soc.AND.PRESENT(sp_alph)) THEN
       CALL soc_sym(&
            &        sym%nop,sym%mrot,sp_beta,sp_alph,cell%amat,&
            &        error)!keep
    ENDIF

    ! determine the group of k
    nopk=0
    ksymloop: DO n=1,sym%nop
       ktest(1)=bk(1)-sym%mrot(1,1,n)*bk(1)-sym%mrot(2,1,n)*bk(2)&
            &        -sym%mrot(3,1,n)*bk(3)
       ktest(2)=bk(2)-sym%mrot(1,2,n)*bk(1)-sym%mrot(2,2,n)*bk(2)&
            &        -sym%mrot(3,2,n)*bk(3)
       ktest(3)=bk(3)-sym%mrot(1,3,n)*bk(1)-sym%mrot(2,3,n)*bk(2)&
            &        -sym%mrot(3,3,n)*bk(3)
       IF (( ABS( ktest(1) - NINT(ktest(1)) ) < eps ) .AND.&
            &         ( ABS( ktest(2) - NINT(ktest(2)) ) < eps ) .AND.&
            &         ( ABS( ktest(3) - NINT(ktest(3)) ) < eps ) .AND.&
            &         (.NOT.error(n)) ) THEN
          nopk=nopk+1
          mrot_k(:,:,nopk)=sym%mrot(:,:,n)
          CYCLE ksymloop
       ENDIF
    ENDDO ksymloop

    DEALLOCATE(error)

    ! Determine the spin-rotations
    ! Double groups not used at the moment, the groups are classified 
    ! without the spin rotations
    !       IF (soc) CALL mrot2su(mrot_k(:,:,1:nopk),amat,su)

    ! identify the group
    ! first determine the classes
    members=1
    nclass=1
    elem(nclass)=1
    belongs(1)=1
    !       IF (soc) THEN 
    ! double group
    !         DO n=1,nopk        
    !            mrot_k(:,:,n+nopk)=mrot_k(:,:,n)
    !            su(:,:,n+nopk)=-su(:,:,n)
    !         ENDDO
    !         nopk=2*nopk
    !       ENDIF

    cloop: DO n=2,nopk
       DO n2=1,nopk
          mtmp=mrot_k(:,:,n2)
          CALL inv3(mrot_k(:,:,n2),mtmpinv,d)
          mtest=MATMUL(MATMUL(mtmp,mrot_k(:,:,n)),mtmpinv)
          !            IF (soc) THEN
          !               sutmp=su(:,:,n2)
          !               sutmpinv=CONJG(TRANSPOSE(sutmp))
          !               sutest=MATMUL(MATMUL(sutmp,su(:,:,n)),sutmpinv)
          !            ENDIF
          DO c=1,nclass
             !               IF (soc) THEN
             !                  IF (ALL((mtest-mrot_k(:,:,elem(c))).EQ.0).AND.
             !c     &   ALL(ABS(REAL(sutest-su(:,:,elem(c)))).LE.0.0001).AND.
             !     &   ALL(ABS(AIMAG(sutest-su(:,:,elem(c)))).LE.0.0001)) THEN
             !                     belongs(n)=c
             !                     CYCLE classloop
             !                  ENDIF
             !               ELSE
             IF (ALL((mtest-mrot_k(:,:,elem(c))).EQ.0)) THEN 
                belongs(n)=c
                members(c)=members(c)+1
                CYCLE cloop
             ENDIF
             !               ENDIF
          ENDDO
       ENDDO
       nclass=nclass+1
       elem(nclass)=n
       belongs(n)=nclass
    ENDDO cloop
    nirr=nclass
    IF (soc) nirr=2*nclass

    ! sort the classes according the rotation type
    mrot2(:,:,1:nopk)=mrot_k(:,:,1:nopk)
    DO c=1,nclass
       mrot_k(:,:,c)=mrot2(:,:,elem(c))
    ENDDO
    !        IF (soc) THEN
    !           su2(:,:,1:nopk)=su(:,:,1:nopk)
    !           DO c=1,nclass
    !              su(:,:,c)=su2(:,:,elem(c))
    !           ENDDO
    !        ENDIF



    ! identify the group
    rot=0
    rot_type=0
    ! rot_type: 2=two fold rot, 3=three fold rot, etc, 7=identity
    !           8=two fold improper rot, ...
    munit=0
    munit(1,1)=1
    munit(2,2)=1
    munit(3,3)=1
    ! determine the number of different rotations
    DO c=1,nclass
       det(c)=&
            &      mrot_k(1,1,c)*mrot_k(2,2,c)*mrot_k(3,3,c)+&
            &      mrot_k(1,2,c)*mrot_k(2,3,c)*mrot_k(3,1,c)+&
            &      mrot_k(2,1,c)*mrot_k(3,2,c)*mrot_k(1,3,c)-&
            &      mrot_k(1,3,c)*mrot_k(2,2,c)*mrot_k(3,1,c)-&
            &      mrot_k(2,3,c)*mrot_k(3,2,c)*mrot_k(1,1,c)-&
            &      mrot_k(2,1,c)*mrot_k(1,2,c)*mrot_k(3,3,c)
       mtest=det(c)*mrot_k(:,:,c)
       rotloop: DO i=1,6
          IF (ALL((mtest-munit).EQ.0)) THEN
             rot(i+(1-det(c))*3)=rot(i+(1-det(c))*3)+1
             IF (i.EQ.1) THEN
                rot_type(c)=7-(1-det(c))*3
             ELSE
                rot_type(c)=i-(1-det(c))*3
             ENDIF
             EXIT rotloop
          ENDIF
          mtest=MATMUL(det(c)*mrot_k(:,:,c),mtest)
       ENDDO rotloop
       IF (ANY((mtest-munit).NE.0)) THEN
          WRITE(6,*) 'grp_k: Cannot find the order of rotation'
       ENDIF
       IF ((rot(5).GT.0).OR.(rot(11).GT.0)) THEN
          WRITE(6,*) 'grp_k: 5 fold rotation found!'
       ENDIF
       CALL euler(c,sym,cell,alpha,beta,gamma)
       CALL rotaxis(alpha,beta,gamma,rax(1:4,c),theta(c))
       IF (soc) THEN
          su(1,1,c)=COS(beta/2.0)*EXP(CMPLX(0.0,-(alpha+gamma)/2.0))
          su(1,2,c)=-SIN(beta/2.0)*EXP(CMPLX(0.0,-(alpha-gamma)/2.0))
          su(2,1,c)=SIN(beta/2.0)*EXP(CMPLX(0.0,(alpha-gamma)/2.0))
          su(2,2,c)=COS(beta/2.0)*EXP(CMPLX(0.0,(alpha+gamma)/2.0))
       ENDIF
    ENDDO

    !<-- Sort the classes
    ! The group elements are sorted in the following way:
    ! First the proper rotations, then improper ones, with increasing rotation angle
    ! Rotations with the same angle are arranged with increasing magnitude of the rotation
    ! axis, e.g. (I, 90 deg around 001, 90 around 110, mirror, ...

    l_sorted=.FALSE.
    DO WHILE (.NOT.l_sorted)
       l_sorted=.TRUE.
       DO c=1,nclass-1
          IF (rot_type(c).LT.rot_type(c+1)) THEN
             mtest=mrot_k(:,:,c)
             mrot_k(:,:,c)=mrot_k(:,:,c+1)
             mrot_k(:,:,c+1)=mtest
             rtmp=theta(c)
             theta(c)=theta(c+1)
             theta(c+1)=rtmp
             rtmp2=rax(:,c)
             rax(:,c)=rax(:,c+1)
             rax(:,c+1)=rtmp2
             itmp=det(c)
             det(c)=det(c+1)
             det(c+1)=itmp
             itmp=rot_type(c)
             rot_type(c)=rot_type(c+1)
             rot_type(c+1)=itmp
             itmp=members(c)
             members(c)=members(c+1)
             members(c+1)=itmp
             IF (soc) THEN
                sutest=su(:,:,c)
                su(:,:,c)=su(:,:,c+1)
                su(:,:,c+1)=sutest
             ENDIF
             l_sorted=.FALSE.
          ENDIF
          IF ((rot_type(c).EQ.rot_type(c+1)).AND.&
               &            (theta(c).GT.theta(c+1))) THEN
             !              IF (theta(c)+(1-det(c))*2.0*pi.GT.
             !     &             theta(c+1)+(1-det(c+1))*2.0*pi) THEN
             mtest=mrot_k(:,:,c)
             mrot_k(:,:,c)=mrot_k(:,:,c+1)
             mrot_k(:,:,c+1)=mtest
             rtmp=theta(c)
             theta(c)=theta(c+1)
             theta(c+1)=rtmp
             rtmp2=rax(:,c)
             rax(:,c)=rax(:,c+1)
             rax(:,c+1)=rtmp2
             itmp=det(c)
             det(c)=det(c+1)
             det(c+1)=itmp
             itmp=members(c)
             members(c)=members(c+1)
             members(c+1)=itmp
             IF (soc) THEN
                sutest=su(:,:,c)
                su(:,:,c)=su(:,:,c+1)
                su(:,:,c+1)=sutest
             ENDIF
             l_sorted=.FALSE.
          ENDIF
          IF ((rot_type(c).EQ.rot_type(c+1)).AND.&
               &            (ABS(theta(c)-theta(c+1)).LT.0.0001).AND.    &
               &            (rax(4,c).GT.rax(4,c+1))) THEN
             !              IF ((ABS(theta(c)+(1-det(c))*2.0*pi-
             !     &           theta(c+1)-(1-det(c+1))*2.0*pi).LT.0.0001).AND.
             !     &           (rax(4,c).GT.rax(4,c+1))) THEN
             mtest=mrot_k(:,:,c)
             mrot_k(:,:,c)=mrot_k(:,:,c+1)
             mrot_k(:,:,c+1)=mtest
             rtmp2=rax(:,c)
             rax(:,c)=rax(:,c+1)
             rax(:,c+1)=rtmp2
             itmp=members(c)
             members(c)=members(c+1)
             members(c+1)=itmp
             IF (soc) THEN
                sutest=su(:,:,c)
                su(:,:,c)=su(:,:,c+1)
                su(:,:,c+1)=sutest
             ENDIF
             l_sorted=.FALSE.
          ENDIF
       ENDDO
    ENDDO
    !>


    WRITE(24,110) bk
    DO c=1,nclass
       IF (det(c).EQ.1) THEN
          WRITE(24,111) NINT(theta(c)*180/pi_const),(rax(1:3,c)),members(c)
       ELSE
          WRITE(24,112) NINT(theta(c)*180/pi_const),(rax(1:3,c)),members(c)
       ENDIF
    ENDDO
110 FORMAT('Symmetry operations of group of k=',3f6.3)
111 FORMAT(i3,1x,'degree proper rotation around ',3f6.2,3x,&
         &      i3,' members in class')
112 FORMAT(i3,1x,'degree improper rotation around ',3f6.2,3x,i3&
         &       ,' members in class')


    !<-- Character tables

    char_table=0.0
    char_table(1,:)=1.0
    grpname='Unknown'
    irrname='Unkno'   ! Only 5 characters wide, cf. ../sympsi.F.
    ! First look the number of classes, within groups with the same number of classes
    ! check the number of different rotations
    SELECT CASE(nclass)
    CASE(1)                           ! only C1
       grpname='C1'
       char_table(1,1)=1.0
       irrname(1)='Gam1'
       IF (soc) THEN
          char_table(2,1)=1.0
          irrname(2)='Gam2'
       ENDIF
    CASE(2)                           ! C-1, C2 and Cs
       IF (rot(2).GT.0) THEN
          grpname='C2'
          char_table(1,1:2)=(/1.0,  1.0/)
          char_table(2,1:2)=(/1.0, -1.0/)
          irrname(1)='Gam1'
          irrname(2)='Gam2'
       ELSE
          grpname='C1h'
          char_table(1,1:2)=(/1.0,  1.0/)
          char_table(2,1:2)=(/1.0, -1.0/)
          irrname(1)='Gam1'
          irrname(2)='Gam2'
          IF (soc) THEN
             char_table(3,1:2)=(/one, -imi/)
             char_table(4,1:2)=(/one,  imi/)
             irrname(3)='Gam3'
             irrname(4)='Gam4'
          ENDIF
       ENDIF
    CASE(3)                           ! C3, D3 and C3v
       IF (ANY(det(1:3).EQ.-1)) THEN
          grpname='C3v'
          char_table(1,1:3)=(/1.0,  1.0,  1.0/)
          char_table(2,1:3)=(/1.0,  1.0, -1.0/)
          char_table(3,1:3)=(/2.0, -1.0,  0.0/)
          irrname(1)='Lam1'
          irrname(2)='Lam2'
          irrname(3)='Lam3'
          IF (soc) THEN
             nirr=6
             char_table(4,1:3)=(/2.0,  1.0,  0.0/)
             char_table(5,1:3)=(/one, -one,  imi/)
             char_table(6,1:3)=(/one, -one, -imi/)
             irrname(4)='Lam6'
             irrname(5)='Lam4'
             irrname(6)='Lam5'
          ENDIF
       ELSE
          grpname='C3'
          omega=EXP(CMPLX(0.0,-2*pi_const/3.0))
          char_table(1,1:3)=(/1.0,  1.0,  1.0/)
          char_table(2,1:3)=(/one, omega, omega**2/)
          char_table(3,1:3)=(/one, omega**2, omega /)
          irrname(1)='Gam1'
          irrname(2)='Gam2'
          irrname(3)='Gam3'
       ENDIF
    CASE(4)                           ! C2h, D2, C2v, C4, S4, and T
       IF((rot(2).EQ.1).AND.(rot(8).EQ.2)) THEN
          grpname='C2v'
          char_table(1,1:4)=(/1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:4)=(/1.0,  1.0, -1.0, -1.0/)
          char_table(3,1:4)=(/1.0, -1.0, -1.0,  1.0/)
          char_table(4,1:4)=(/1.0, -1.0,  1.0, -1.0/)
          irrname(1)='Z1'
          irrname(2)='Z2'
          irrname(3)='Z3'
          irrname(4)='Z4'            
       ELSE IF (rot(3).GT.1) THEN
          grpname='T'
          omega=EXP(CMPLX(0.0,-2*pi_const/3.0))
          char_table(1,1:4)=(/1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:4)=(/one,  one, omega, omega**2/)
          char_table(3,1:4)=(/one,  one, omega**2,  omega/)
          char_table(4,1:4)=(/3.0, -1.0,  0.0,  0.0/)
          irrname(1)='Gam1'
          irrname(2)='Gam2'
          irrname(3)='Gam3'
          irrname(4)='Gam4'
       ELSE IF ((rot(4).GT.1).OR.(rot(10).GT.1)) THEN
          IF (rot(4).GT.1) grpname='C4'
          IF (rot(10).GT.1) grpname='S4'
          char_table(1,1:4)=(/1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:4)=(/1.0, -1.0,  1.0, -1.0/)
          char_table(3,1:4)=(/one, -imi, -one,  imi/)
          char_table(4,1:4)=(/one,  imi, -one, -imi/)
          irrname(1)='Gam1'
          irrname(2)='Gam2'
          irrname(3)='Gam3'
          irrname(4)='Gam4'
          IF (soc) THEN
             nirr=8
             omega=EXP(CMPLX(0.0,-pi_const/4.0))
             char_table(5,1:4)=(/one, omega, -imi, &
                  &                   -CONJG(omega)/)
             char_table(6,1:4)=(/one, CONJG(omega), imi, &
                  &                   -omega/)
             char_table(7,1:4)=(/one, -omega, -imi, &
                  &                   CONJG(omega)/)
             char_table(8,1:4)=(/one, -CONJG(omega), imi, &
                  &                   omega/)
             irrname(5)='Gam5'
             irrname(6)='Gam6'
             irrname(7)='Gam7'
             irrname(8)='Gam8'                   
          ENDIF
       ELSE IF ((rot(2).EQ.1).AND.(rot(8).EQ.1)) THEN
          grpname='ZKUS'
       ENDIF
    CASE(5)                           ! D4, C4v, D2d, O and Td
       IF (rot(3).GT.0) THEN
          grpname='Td'
          char_table(1,1:5)=(/1.0,  1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:5)=(/1.0, -1.0,  1.0, -1.0,  1.0/)
          char_table(3,1:5)=(/2.0,  0.0,  2.0,  0.0, -1.0/)
          char_table(4,1:5)=(/3.0,  1.0, -1.0, -1.0,  0.0/)
          char_table(5,1:5)=(/3.0, -1.0, -1.0,  1.0,  0.0/)
          irrname(1)='P1'
          irrname(2)='P2'
          irrname(3)='P3'
          irrname(4)='P5'
          irrname(5)='P4'
       ELSE IF (rot(10).GT.0) THEN
          grpname='D2d'
          char_table(1,1:5)=(/1.0,  1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:5)=(/1.0,  1.0, -1.0,  1.0, -1.0/)
          char_table(3,1:5)=(/1.0,  1.0,  1.0, -1.0, -1.0/)
          char_table(4,1:5)=(/1.0,  1.0, -1.0, -1.0,  1.0/)
          char_table(5,1:5)=(/2.0, -2.0,  0.0,  0.0,  0.0/)
          !               standard order: E, 2, 2_x, -4, m
          irrname(1)='W1'
          irrname(2)='W2'
          irrname(3)='W1`'
          irrname(4)='W2`'
          irrname(5)='W3'
       ELSE
          grpname='C4v'
          char_table(1,1:5)=(/1.0,  1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:5)=(/1.0,  1.0,  1.0, -1.0, -1.0/)
          char_table(3,1:5)=(/1.0, -1.0,  1.0,  1.0, -1.0/)
          char_table(4,1:5)=(/1.0, -1.0,  1.0, -1.0,  1.0/)
          char_table(5,1:5)=(/2.0,  0.0, -2.0,  0.0,  0.0/)
          irrname(1)='Del1'
          irrname(2)='Del1`'
          irrname(3)='Del2'
          irrname(4)='Del2`'
          irrname(5)='Del5'
       ENDIF
    CASE(6)                           ! C3i, D3d, C6, C3h, D6, C6v and D3h
       IF (rot(2).EQ.0) THEN
          IF (rot(7).EQ.0) THEN
             grpname='C3h'
             char_table(:,:) = 0.0
             WRITE(24,*) 'C3h: Character table missing in grp_k.F'
          ELSE
             grpname='C3i'
             omega=CMPLX(-0.5,SQRT(3./4.))
             char_table(1,1:6)=(/one, one, one, one, one, one/)
             char_table(2,1:6)=(/one, one, one,-one,-one,-one/)
             char_table(3,1:6)=(/one, omega, CONJG(omega),&
                  &                                one, omega, CONJG(omega)/)
             char_table(4,1:6)=(/one, omega, CONJG(omega),&
                  &                               -one,-omega,-CONJG(omega)/)
             char_table(5,1:6)=(/one, CONJG(omega), omega,&
                  &                                one, CONJG(omega), omega/)
             char_table(6,1:6)=(/one, CONJG(omega), omega,&
                  &                               -one,-CONJG(omega),-omega/)
             irrname(1)='Ag  '
             irrname(2)='Au  '
             irrname(3)='E1g '
             irrname(4)='E1u '
             irrname(5)='E2g '
             irrname(6)='E2u '
             IF (soc) THEN
                char_table(7,1:6)=(/one,-one, one, one,-one, one/)
                char_table(8,1:6)=(/one,-one, one,-one, one,-one/)
                char_table( 9,1:6)=(/one,-omega, CONJG(omega),&
                     &                                one,-omega, CONJG(omega)/)
                char_table(10,1:6)=(/one,-omega, CONJG(omega),&
                     &                               -one, omega,-CONJG(omega)/)
                char_table(11,1:6)=(/one,-CONJG(omega), omega,&
                     &                                one,-CONJG(omega), omega/)
                char_table(12,1:6)=(/one,-CONJG(omega), omega,&
                     &                               -one, CONJG(omega),-omega/)
                irrname(7)="Ag' "
                irrname(8)="Au' "
                irrname(9)="E1g'"
                irrname(10)="E1u'"
                irrname(11)="E2g'"
                irrname(12)="E2u'"
             ENDIF
          ENDIF
       ELSEIF (rot(2).GT.1) THEN
          grpname='D6'
          char_table(1,1:6)=(/1.0,  1.0,  1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:6)=(/1.0,  1.0,  1.0,  1.0, -1.0, -1.0/)
          char_table(3,1:6)=(/1.0, -1.0,  1.0, -1.0,  1.0, -1.0/)
          char_table(4,1:6)=(/1.0, -1.0,  1.0, -1.0, -1.0,  1.0/)
          char_table(5,1:6)=(/2.0,  1.0, -1.0, -2.0,  0.0,  0.0/)
          char_table(6,1:6)=(/2.0, -1.0, -1.0,  2.0,  0.0,  0.0/)
          irrname(1)='Gam1'
          irrname(2)='Gam2'
          irrname(3)='Gam3'
          irrname(4)='Gam4'
          irrname(5)='Gam6'
          irrname(6)='Gam5'
       ELSE
          IF (rot(6).EQ.0) THEN
             grpname='D3d'
             char_table(1,1:6)=(/1.0,  1.0,  1.0,  1.0,  1.0,  1.0/)
             char_table(2,1:6)=(/1.0,  1.0, -1.0,  1.0,  1.0, -1.0/)
             char_table(3,1:6)=(/2.0, -1.0,  0.0,  2.0, -1.0,  0.0/)
             char_table(4,1:6)=(/1.0,  1.0,  1.0, -1.0, -1.0, -1.0/)
             char_table(5,1:6)=(/1.0,  1.0, -1.0, -1.0, -1.0,  1.0/)
             char_table(6,1:6)=(/2.0, -1.0,  0.0, -2.0,  1.0,  0.0/)
             irrname(1)='A1g'
             irrname(2)='A2g'
             irrname(3)='Eg'
             irrname(4)='A1u'
             irrname(5)='A2u'
             irrname(6)='Eu'
          ELSEIF (rot(6).GT.1) THEN
             grpname='C6' 
             char_table(:,:) = 0.0
             WRITE(24,*) 'C6: Character table missing in grp_k.F'
          ELSE
             IF (members(3).GT.1) grpname='D3h'    ! maybe this works
             IF (members(3).EQ.1) grpname='C6v'
             char_table(1,1:6)=(/1.0,  1.0,  1.0,  1.0,  1.0,  1.0/)
             char_table(2,1:6)=(/1.0,  1.0,  1.0,  1.0, -1.0, -1.0/)
             char_table(3,1:6)=(/1.0, -1.0,  1.0, -1.0,  1.0, -1.0/)
             char_table(4,1:6)=(/1.0, -1.0,  1.0, -1.0, -1.0,  1.0/)
             char_table(5,1:6)=(/2.0,  1.0, -1.0, -2.0,  0.0,  0.0/)
             char_table(6,1:6)=(/2.0, -1.0, -1.0,  2.0,  0.0,  0.0/)
             irrname(1)='Gam1'
             irrname(2)='Gam2'
             irrname(3)='Gam3'
             irrname(4)='Gam4'
             irrname(5)='Gam6'
             irrname(6)='Gam5'
          ENDIF
       ENDIF
    CASE(8)                           ! Th, C4h and D2h
       IF (rot(3).GT.0) THEN
          grpname='Th'
          char_table(:,:) = 0.0
          WRITE(24,*) 'Th: Character table missing in grp_k.F'
       ELSE IF (rot(2).EQ.3) THEN
          grpname='D2h'
          char_table(:,:) = 0.0
          WRITE(24,*) 'D2h: Character table missing in grp_k.F'
       ELSE
          grpname='C4h'
          char_table(:,:) = 0.0
          WRITE(24,*) 'C4h: Character table missing in grp_k.F'
       ENDIF
    CASE(10)                           ! D4h and Oh
       IF (rot(3).GT.0) THEN
          grpname='Oh'
          char_table(1,1:10)= (/1.0,  1.0,  1.0,  1.0,  1.0,&
               &                                1.0,  1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:10)= (/1.0, -1.0,  1.0,  1.0, -1.0,&
               &                                1.0, -1.0,  1.0,  1.0, -1.0/)
          char_table(3,1:10)= (/2.0,  0.0, -1.0,  2.0,  0.0,&
               &                                2.0,  0.0, -1.0,  2.0,  0.0/)
          char_table(4,1:10)= (/3.0,  1.0,  0.0, -1.0, -1.0,&
               &                                3.0,  1.0,  0.0, -1.0, -1.0/)
          char_table(5,1:10)= (/3.0, -1.0,  0.0, -1.0,  1.0,&
               &                                3.0, -1.0,  0.0, -1.0,  1.0/)
          char_table(6,1:10)= (/1.0,  1.0,  1.0,  1.0,  1.0,&
               &                               -1.0, -1.0, -1.0, -1.0, -1.0/)
          char_table(7,1:10)= (/1.0, -1.0,  1.0,  1.0, -1.0,&
               &                               -1.0,  1.0, -1.0, -1.0,  1.0/)
          char_table(8,1:10)= (/2.0,  0.0, -1.0,  2.0,  0.0,&
               &                               -2.0,  0.0,  1.0, -2.0,  0.0/)
          char_table(9,1:10)= (/3.0,  1.0,  0.0, -1.0, -1.0,&
               &                               -3.0, -1.0,  0.0,  1.0,  1.0/)
          char_table(10,1:10)=(/3.0, -1.0,  0.0, -1.0,  1.0,&
               &                               -3.0,  1.0,  0.0,  1.0, -1.0/)

          irrname(1)='Gam1+'
          irrname(2)='Gam2+'
          irrname(3)='Gam3+'
          irrname(4)='Gam4+'
          irrname(5)='Gam5+'
          irrname(6)='Gam1-'
          irrname(7)='Gam2-'
          irrname(8)='Gam3-'
          irrname(9)='Gam4-'
          irrname(10)='Gam5-'
       ELSE
          grpname='D4h'
          char_table(1,1:10)= (/1.0,  1.0,  1.0,  1.0,  1.0,&
               &                                1.0,  1.0,  1.0,  1.0,  1.0/)
          char_table(2,1:10)= (/1.0,  1.0,  1.0, -1.0, -1.0,&
               &                                1.0,  1.0,  1.0, -1.0, -1.0/)
          char_table(3,1:10)= (/1.0, -1.0,  1.0,  1.0, -1.0,&
               &                                1.0, -1.0,  1.0,  1.0, -1.0/)
          char_table(4,1:10)= (/1.0, -1.0,  1.0, -1.0,  1.0,&
               &                                1.0, -1.0,  1.0, -1.0,  1.0/)
          char_table(5,1:10)= (/2.0,  0.0, -2.0,  0.0,  0.0,&
               &                                2.0,  0.0, -2.0,  0.0,  0.0/)
          char_table(6,1:10)= (/1.0,  1.0,  1.0,  1.0,  1.0,&
               &                               -1.0, -1.0, -1.0, -1.0, -1.0/)
          char_table(7,1:10)= (/1.0,  1.0,  1.0, -1.0, -1.0,&
               &                               -1.0, -1.0, -1.0,  1.0,  1.0/)
          char_table(8,1:10)= (/1.0, -1.0,  1.0,  1.0, -1.0,&
               &                               -1.0,  1.0, -1.0, -1.0,  1.0/)
          char_table(9,1:10)= (/1.0, -1.0,  1.0, -1.0,  1.0,&
               &                               -1.0,  1.0, -1.0,  1.0, -1.0/)
          char_table(10,1:10)=(/2.0,  0.0, -2.0,  0.0,  0.0,&
               &                               -2.0,  0.0,  2.0,  0.0,  0.0/)

          ! standard char-table assumes the order: E, 4, 2, 2_h, 2_h' where 4 and 2 
          ! rotate around the same axis and 2_h is perpendicular to this axis. Check:

          CALL cross( rax(1,2),rax(1,3),kt )
          IF (kt(1)**2+kt(2)**2+kt(3)**2 < eps) THEN
             !                  all ok
          ELSE 
             CALL cross( rax(1,2),rax(1,4),kt )
             IF (kt(1)**2+kt(2)**2+kt(3)**2 < eps) THEN      ! change 3 & 4
                CALL change_column(char_table,10,10,3,4)
                ! check, whether element 3 is perpendicular to element 2
                IF (DOT_PRODUCT(rax(1:3,2),rax(1:3,3)) < eps) THEN
                   !                     all ok
                ELSE             ! change 4 & 5
                   WRITE(*,*) DOT_PRODUCT(rax(1:3,2),rax(1:3,3))
                   CALL change_column(char_table,10,10,4,5)
                ENDIF
             ELSE
                CALL cross( rax(1,2),rax(1,5),kt )
                IF (kt(1)**2+kt(2)**2+kt(3)**2 < eps) THEN    ! change 3 & 5
                   CALL change_column(char_table,10,10,3,5)
                   IF (DOT_PRODUCT(rax(1:3,2),rax(1:3,3)) < eps) THEN ! change 4 & 5
                      CALL change_column(char_table,10,10,4,5)
                   ELSE             
                      !                     all ok
                   ENDIF
                ELSE
                   CALL juDFT_error("D4h",calledby="grp_k")
                ENDIF
             ENDIF
          ENDIF
          irrname(1)='A1g'
          irrname(2)='A2g'
          irrname(3)='B1g'
          irrname(4)='B2g'
          irrname(5)='Eg'
          irrname(6)='A1u'
          irrname(7)='A2u'
          irrname(8)='B1u'
          irrname(9)='B2u'
          irrname(10)='Eu'
       ENDIF
    CASE(12)                           ! D6h and C6h
       IF (rot(2).EQ.3) THEN
          grpname='D6h'
          char_table(:,:) = 0.0
          WRITE(24,*) 'D6h: Character table missing in grp_k.F'
       ELSE
          grpname='C6h'
          char_table(:,:) = 0.0
          WRITE(24,*) 'C6h: Character table missing in grp_k.F'
       ENDIF
    CASE DEFAULT
       WRITE(24,*) 'Group of k not identified'
    END SELECT

    !>

  END SUBROUTINE grp_k
  !-------------------------------------------------------------------------------

  !<--SUBROUTINE rotaxis(alpha,beta,gamma,rax,theta)
  SUBROUTINE rotaxis(alpha,beta,gamma,rax,theta)

    ! Determines the rotation axis based on the Euler angles
    IMPLICIT NONE
    REAL, INTENT(IN) :: alpha,beta,gamma
    REAL, INTENT(OUT) :: rax(4),theta 

    REAL :: costhe2,the2
    INTEGER :: i

    costhe2=COS(beta/2.0)*COS((alpha+gamma)/2.0)
    the2=ACOS(costhe2)
    rax=0.0
    IF (the2.GT.0.00001) THEN
       rax(1)=-SIN(beta/2.0)*SIN((alpha-gamma)/2.0)/SIN(the2)
       rax(2)= SIN(beta/2.0)*COS((alpha-gamma)/2.0)/SIN(the2)
       rax(3)= COS(beta/2.0)*SIN((alpha+gamma)/2.0)/SIN(the2)
    ENDIF
    DO i=1,3
       IF (ABS(rax(i)).GT.0.0001) THEN
          rax(i)=rax(i)/ABS(rax(i))
       ENDIF
    ENDDO
    rax(4)=rax(1)**2+rax(2)**2+rax(3)**2
    theta=the2*2.0
  END SUBROUTINE rotaxis
  !>

  !<--Subroutine euler
  SUBROUTINE euler(n,sym,cell,alpha,beta,gamma)

    ! determines the Euler angles corresponding the proper rotation part of mrot      

    USE m_constants, ONLY : pi_const
    USE m_inv3
    USE m_types
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: n
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_cell),INTENT(IN)    :: cell
    
    REAL, INTENT(OUT) :: alpha,beta,gamma

    INTEGER :: det
    REAL :: mprop(3,3)      
    REAL :: sina,sinb,sinc,cosa,cosb,cosc
    REAL :: amatinv(3,3),detr

    CALL inv3(cell%amat,amatinv,detr)
    det= sym%mrot(1,1,n)*sym%mrot(2,2,n)*sym%mrot(3,3,n) +&
         &        sym%mrot(1,2,n)*sym%mrot(2,3,n)*sym%mrot(3,1,n) +&
         &        sym%mrot(2,1,n)*sym%mrot(3,2,n)*sym%mrot(1,3,n) -&
         &        sym%mrot(1,3,n)*sym%mrot(2,2,n)*sym%mrot(3,1,n) -&
         &        sym%mrot(2,3,n)*sym%mrot(3,2,n)*sym%mrot(1,1,n) -&
         &        sym%mrot(2,1,n)*sym%mrot(1,2,n)*sym%mrot(3,3,n)

    ! Take the proper rotation         
    mprop=REAL(det)*MATMUL(cell%amat,MATMUL(REAL(sym%mrot(:,:,n)),amatinv))
    ! Euler angles         
    cosb = mprop(3,3)
    sinb = 1.00 - cosb*cosb
    sinb = MAX(sinb,0.00)
    sinb = SQRT(sinb)
    !
    ! if beta = 0 or pi , only alpha+gamma or -gamma have a meaning:
    !
    IF ( ABS(sinb).LT.1.0e-5 ) THEN 
       beta = 0.0
       IF ( cosb.LT.0.0 ) beta = pi_const
       gamma = 0.0
       cosa = mprop(1,1)/cosb
       sina = mprop(1,2)/cosb
       IF ( ABS(sina).LT.1.0e-5 ) THEN
          alpha=0.0
          IF ( cosa.LT.0.0 ) alpha=alpha+pi_const
       ELSE
          alpha = 0.5*pi_const - ATAN(cosa/sina)
          IF ( sina.LT.0.0 ) alpha=alpha+pi_const
       ENDIF
    ELSE
       beta = 0.5*pi_const - ATAN(cosb/sinb)
       !     
       ! determine alpha and gamma from d13 d23 d32 d31
       !
       cosa = mprop(3,1)/sinb
       sina = mprop(3,2)/sinb
       cosc =-mprop(1,3)/sinb
       sinc = mprop(2,3)/sinb
       IF ( ABS(sina).LT.1.0e-5 ) THEN
          alpha=0.0
          IF ( cosa.LT.0.0 ) alpha=alpha+pi_const
       ELSE
          alpha = 0.5*pi_const - ATAN(cosa/sina)
          IF ( sina.LT.0.0 ) alpha=alpha+pi_const
       ENDIF
       IF ( ABS(sinc).LT.1.0e-5 ) THEN
          gamma = 0.0
          IF ( cosc.LT.0.0 ) gamma=gamma+pi_const
       ELSE
          gamma = 0.5*pi_const - ATAN(cosc/sinc)
          IF ( sinc.LT.0.0 ) gamma=gamma+pi_const
       ENDIF

    ENDIF

  END SUBROUTINE euler
  !>
  !-----------------------------------------------------------
  SUBROUTINE change_column( char_table,nx,ny,r1,r2 )

    INTEGER, INTENT (IN)   :: nx,ny,r1,r2
    COMPLEX, INTENT(INOUT) :: char_table(:,:)

    COMPLEX tmpc(nx)

    tmpc(1:nx) = char_table(1:nx,r1)
    char_table(1:nx,r1) = char_table(1:nx,r2)
    char_table(1:nx,r2) = tmpc(1:nx)
    tmpc(1:nx) = char_table(1:nx,r1+nx/2)
    char_table(1:nx,r1+nx/2) = char_table(1:nx,r2+nx/2)
    char_table(1:nx,r2+nx/2) = tmpc(1:nx)

  END SUBROUTINE change_column
END MODULE m_grp_k
