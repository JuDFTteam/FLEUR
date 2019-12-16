MODULE m_symproperties
  USE m_juDFT
  !********************************************************************
  !      calculates various properties about each symmetry operation
  !      and the space group + lattice
  !********************************************************************
CONTAINS
  SUBROUTINE symproperties(&
       optype,oldfleur,nops,multtab,amat,&
       symor,mrot,tau,invsym,invs,zrfs,invs2,nop,nop2)

    IMPLICIT NONE

    !===> Arguments
    INTEGER, INTENT (IN)  :: nops
    INTEGER, INTENT (IN)  :: optype(:),multtab(nops,nops)
    REAL,    INTENT (IN)  :: amat(3,3)
    LOGICAL, INTENT (IN)  :: oldfleur
    LOGICAL, INTENT (OUT) :: invsym,invs,zrfs,invs2
    INTEGER, INTENT (OUT) :: nop,nop2                ! if .oldfleur. nop <=nops
    LOGICAL, INTENT (INOUT) :: symor
    INTEGER, INTENT (INOUT) :: mrot(3,3,nops)
    REAL,    INTENT (INOUT) :: tau(3,nops)

    !===> Local Variables

    INTEGER invsop, zrfsop, invs2op, magicinv
    INTEGER i,j,na,nn,n,dbgfh
    INTEGER indtwo(SIZE(optype)), usedop(SIZE(optype))
    INTEGER mrotaux(3,3,SIZE(optype))
    INTEGER mtab(nops,nops),iop(nops)
    REAL    tauaux(3,SIZE(optype)), eps12
    LOGICAL zorth           ! true, if z-axis is othorgonal

    invsym  = .FALSE.
    invs    = .FALSE.
    zrfs    = .FALSE.
    invs2   = .FALSE.
    invsop  = 0
    zrfsop  = 0
    invs2op = 0
    eps12   = 1.0e-12
    dbgfh   = 6

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    symor= .TRUE.
    IF ( ANY ( ABS(tau(:,1:nops)) > eps12 ) ) symor = .FALSE.

    DO n = 1, nops

       IF ( optype(n) == -1 )THEN
          invsym = .TRUE.
          !--->       check if we have inversion as a symmetry operation
          IF ( ALL ( ABS( tau(:,n) ) < eps12 ) ) THEN
             invsop = n
             invs = .TRUE.
          ENDIF
       ENDIF

       IF ( optype(n) == -2 )THEN
          !--->      check for z-reflection
          IF ( mrot(3,3,n) == -1 .AND. ALL(ABS(tau(:,n))<eps12) ) THEN
             zrfsop = n
             zrfs = .TRUE.
          ENDIF
       ENDIF

       IF ( optype(n) == 2 )THEN
          !--->      check for 2d inversion
          IF ( mrot(3,3,n) == 1 .AND. ALL(ABS(tau(:,n))<eps12) ) THEN
             invs2op = n
             invs2 = .TRUE.
          ENDIF
       ENDIF

    ENDDO !nops

    IF ( amat(3,1)==0.00 .AND. amat(3,2)==0.00 .AND.amat(1,3)==0.00 .AND. amat(2,3)==0.00 ) THEN
       zorth= .TRUE.
    ELSE       
       zorth= .FALSE.
       ! reset the following...
       zrfs    = .FALSE.
       invs2   = .FALSE.
    ENDIF

    WRITE(6,*)
    WRITE(6,*) 'DBG: symor,zorth,oldfleur :', symor,zorth,oldfleur
    WRITE(6,'(1x,a13,48i5)') 'DBG: optype :', optype(1:nops)
    WRITE(6,*) 'DBG: invsym,invs,zrfs,invs2 :', invsym,invs,zrfs,invs2
    WRITE(6,'(1x,a45,3i5)') 'DBG: (before reorder) invsop,zrfsop,invs2op :', invsop,zrfsop,invs2op

    IF ( (.NOT.oldfleur) .OR. (.NOT.zorth) ) THEN
       nop = nops
       nop2 = 0
       IF (.NOT.oldfleur) RETURN
    ENDIF
    IF ( oldfleur .AND. (.NOT.zorth) ) THEN
       CALL juDFT_error("oldfleur = t and z-axis not orthogonal",calledby ="symproperties")
    ENDIF
    nop = nops

    !---> now we have to sort the ops to find the two-dimensional ops
    !---> and their 3-dim inverted or z-reflected counterparts

    mrotaux(:,:,1:nops) = mrot(:,:,1:nops)
    tauaux(:,1:nops) = tau(:,1:nops)

    DO i=1,nops
       indtwo(i)= i
    ENDDO

    nop2=0
    DO i = 1, nops
       IF ( mrot(3,3,i) == 1 ) THEN
          nop2 = nop2 + 1
          indtwo(nop2)= i
       ENDIF
    ENDDO

    !dbg  write(dbgfh,*) 'DBG: nop2 : ', nop2

    magicinv = 0
    IF (zrfs) magicinv = zrfsop
    IF (invs) magicinv = invsop
    usedop = 1

    IF ( magicinv > 0 ) THEN
       DO i = 1, nop2
          j = indtwo(i)
          mrot(:,:,i) = mrotaux(:,:,j)
          tau(:,i)    = tauaux(:,j)
          usedop(j) = usedop(j) - 1
          j = multtab(magicinv,indtwo(i))
          mrot(:,:,i+nop2) =  mrotaux(:,:,j)
          tau(:,i+nop2) = tauaux(:,j)
          usedop(j) = usedop(j) - 1
       ENDDO

       IF ( ANY( usedop(1:nops) < 0 ) ) THEN
          WRITE (dbgfh,*) 'DBG: usedop : ', usedop(1:nops)
          CALL juDFT_error("Fatal Error! #01",calledby="symproperties")
       ENDIF

       nop = 2*nop2
       IF ( nop.NE.nops ) THEN
          n = 0
          DO i = 1, nops
             IF ( usedop(i) == 1 ) THEN
                n = n + 1
                mrot(:,:,nop+n) =  mrotaux(:,:,i)
                tau(:,nop+n) = tauaux(:,i)
             ENDIF
          ENDDO

          !dbg      write(dbgfh,*) 'DBG: nops, nop, n : ', nops, nop, n

          IF ( n+nop /= nops )  CALL juDFT_error("Fatal Error! #02" ,calledby ="symproperties")
       ENDIF

    ENDIF


    !---> check for nonsymmorphic translations in z-direction in
    !---> the film (oldfleur=t) case

    IF ( oldfleur ) THEN

       n = 1
       DO WHILE (n <= nop)
          IF (ABS(tau(3,n)) > 0.000001) THEN
             mrotaux(:,:,1) = mrot(:,:,n)
             tauaux(:,1) = tau(:,n)
             DO nn = n+1, nops
                mrot(:,:,nn-1) = mrot(:,:,nn)
                tau(:,nn-1) = tau(:,nn)
             ENDDO
             mrot(:,:,nops) = mrotaux(:,:,1)
             tau(:,nops) = tauaux(:,1) 
             nop = nop - 1
             WRITE(*,*) 'op',n,'removed'
          ELSE
             n = n + 1
          ENDIF
       ENDDO
       WRITE(*,*) 'nop =',nop

    ENDIF

    IF ( oldfleur .AND. nop.NE.nops ) THEN
       WRITE(6,'(/," Full space group has",i3," operations.",/)') nops
       WRITE(6,'(i3," operations violate the 2d symmetry in fleur21"," and have been removed.",/)') nops-nop
       DO n = nop+1, nops
          WRITE(6,'(" operation",i3,":  ")') n
          WRITE(6,'(15x,"  (",3i3," )  (",f6.3," )")') ((mrot(j,i,n),i=1,3),tau(j,n),j=1,3)
       ENDDO
       WRITE(6,'(/,"Reduced space group has",i3," operations.",/)') nop
       !        nops = nop
    ELSE
       nop = nops
    ENDIF

  END SUBROUTINE symproperties
END MODULE m_symproperties
