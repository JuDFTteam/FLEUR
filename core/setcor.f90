MODULE m_setcor
  USE m_juDFT
CONTAINS
  SUBROUTINE setcor(itype,jspins,atoms,bmu,nst,kappa,nprnc,occ)
    !
    !     *****************************************************
    !     sets the values of kappa and occupation numbers of
    !     the neutral atoms.
    !         following code by m. weinert  february 1982
    !     *****************************************************

    USE m_types
    USE m_types
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !
    !     .. Scalar Arguments ..
    INTEGER,INTENT (IN) :: itype,jspins
    REAL,INTENT (INOUT) :: bmu
    !     ..
    INTEGER,INTENT (OUT) :: nst
    INTEGER,INTENT (OUT) :: kappa(:),nprnc(:)
    REAL,INTENT (OUT)    :: occ(:,:)
    !     ..
    !     .. Local Scalars ..
    INTEGER iz,jz,jz0,k,n,jspin
    INTEGER k_h(2),n_h(2)
    REAL fj,l,bmu_l,o_h(2), fac(2)
    LOGICAL l_clf
    CHARACTER(len=13) :: fname
    !     ..

    l_clf = .FALSE.  
    WRITE(fname,"('corelevels.',i2.2)") NINT(atoms%zatom(itype))
    INQUIRE (file=fname, exist=l_clf)

    IF (l_clf) THEN
       OPEN (61,file=fname,form='formatted')
       READ (61,'(i3)') nst
       IF (bmu.LT.0.001) bmu = 999.
       IF (nst > size(kappa))  CALL juDFT_error("corelevels: nst > nstd" ,calledby ="setcor")
       DO n = 1, nst
          fac(1) = 1.0 ; fac(2) = 1.0
          READ (61,'(4i3)') nprnc(n),kappa(n),n_h(1),n_h(2)
          IF ( n_h(1) < 0 )  fac(1) = -0.5
          IF ( n_h(2) < 0 )  fac(2) = -0.5
          IF (jspins.EQ.1) THEN
             occ(n,1) = fac(1) * n_h(1) + fac(2) * n_h(2) 
          ELSE
             occ(n,1) = fac(1) * n_h(1) ; occ(n,2) = fac(2) * n_h(2)
          ENDIF
          !         write(*,*) nprnc(n),kappa(n),occ(n,1), occ(n,2)
       ENDDO
       CLOSE (61)
       RETURN
    ELSE
       jspin=1
    ENDIF

    IF (atoms%zatom(itype)>92.01e0)  CALL juDFT_error(" z > 92",calledby ="setcor"&
         &     )

    jz0 = atoms%zatom(itype) + 0.01e0
    jz = jz0
    k = 0
    DO n = 1,7 !maximal main quantuum number =7
       IF (jz.LE.0) EXIT
       !--->    s states
       k = k + 1
       nprnc(k) = n
       kappa(k) = -1
       occ(k,1) = 2
       jz = jz - 2
       IF (jz.LE.0) EXIT
       !--->    p states
       IF (n.EQ.1) CYCLE
       k = k + 1
       nprnc(k) = n
       kappa(k) = 1
       occ(k,1) = 2
       jz = jz - 2
       IF (jz.LE.0) EXIT
       k = k + 1
       nprnc(k) = n
       kappa(k) = -2
       occ(k,1) = 4
       jz = jz - 4
       IF (jz.LE.0) EXIT
       !--->    d functions
       iz = 0
       IF (n.EQ.3 .AND. jz0.GT.20) iz = MIN(jz0-20,4)
       IF (n.EQ.4 .AND. jz0.GT.38) iz = MIN(jz0-38,4)
       IF (n.EQ.4 .AND. jz0.EQ.41) iz = 4
       IF (n.EQ.5 .AND. jz0.GT.70) iz = MIN(jz0-70,4)
       IF (n.EQ.5 .AND. (jz0.EQ.57.OR.jz0.EQ.64)) iz = 1
       IF (n.EQ.6 .AND. jz0.GT.88) iz = 1
       IF (n.EQ.6 .AND. jz0.EQ.90) iz = 2
       IF (iz.NE.0) THEN 
          k = k + 1
          nprnc(k) = n
          kappa(k) = 2
          occ(k,1) = iz
          jz = jz - iz
          IF ((n==6).AND.(iz.GE.4)) CYCLE
          IF (iz.GE.4 .AND. .NOT.(n.EQ.4 .AND. jz0.EQ.41) .AND. .NOT. (n.EQ.5 .AND. jz0.EQ.74)) THEN
             iz = 1
             IF (n.EQ.3 .AND. jz0.GT.25) iz = MIN(jz0-24,6)
             IF (n.EQ.3 .AND. jz0.EQ.29) iz = 6
             IF (n.EQ.4 .AND. jz0.GT.43) iz = jz0 - 41
             IF (n.EQ.4 .AND. jz0.GT.45) iz = 6
             IF (n.EQ.5 .AND. jz0.GT.75) iz = jz0 - 74
             IF (n.EQ.5 .AND. jz0.GT.77) iz = 6
             k = k + 1
             nprnc(k) = n
             kappa(k) = -3
             occ(k,1) = iz
             jz = jz - iz
             
          ENDIF
       ENDIF
       !--->    f states
       IF (n==4) THEN
          !+gu  IF (jz0.LE.57) GO TO 50
          IF (jz0.LE.62) CYCLE
          k = k + 1
          nprnc(k) = n
          kappa(k) = 3
          iz = 6
          IF (jz0.LT.62) THEN
             iz = jz0 - 56
             occ(k,1) = iz
             jz = jz - iz
             CYCLE
          ENDIF
          occ(k,1) = iz
          jz = jz - iz
          iz = 8
          k = k + 1
          nprnc(k) = n
          kappa(k) = -4
          IF (jz0.LT.70) THEN
             iz = jz0 - 62
             IF (jz0.EQ.64) iz = 1
             occ(k,1) = iz
             jz = jz - iz
             CYCLE
          ENDIF
          occ(k,1) = iz
          jz = jz - iz
       ENDIF
       IF (n.NE.5) CYCLE
       IF (jz0.LE.90) CYCLE
       k = k + 1
       nprnc(k) = n
       kappa(k) = 3
       iz = jz0 - 89
       occ(k,1) = iz
       jz = jz - iz
    ENDDO
    nst = k
    IF (k.GE.1) occ(k,1) = occ(k,1) + jz
    !
    ! add magnetic moments
    !
    IF (jspins.EQ.2) THEN
       bmu_l = bmu
       DO k = 1,nst
          occ(k,jspins) = occ(k,1)/2.0 
          occ(k,1) = occ(k,jspins)  
       ENDDO
       kloop:DO k = nst,1,-1
          fj = iabs(kappa(k)) - 0.5e0
          l = fj + 0.5e0*isign(1,kappa(k)) + 0.01e0
          ! polarize (d,f) only
          IF (l.GT.1.99) THEN
             IF (2*occ(k,1).GE.ABS(bmu_l)) THEN
                occ(k,1) = occ(k,1) + bmu_l/2.
                occ(k,jspins) = occ(k,jspins) - bmu_l/2.
                EXIT kloop
             ELSE
                IF (bmu_l.GT.0) THEN
                   occ(k,1) = 2.0*occ(k,1)
                   occ(k,jspins) = 0.0
                   bmu_l = bmu_l - occ(k,1)
                ELSE
                   occ(k,jspins) = 2.0*occ(k,jspins)
                   occ(k,1) = 0.0
                   bmu_l = bmu_l + occ(k,jspins)
                ENDIF
             ENDIF
          ENDIF
       ENDDO kloop
    ENDIF
    !
    IF (atoms%zatom(itype).EQ.65) THEN
       k_h(1) = kappa(15) ; n_h(1) = nprnc(15) ; o_h(1) = occ(15,1)
       k_h(2) = kappa(16) ; n_h(2) = nprnc(16) ; o_h(2) = occ(16,1)
       kappa(15)= kappa(17) ; nprnc(15)=nprnc(17) ; occ(15,1)=occ(17,1)
       kappa(16)= kappa(18) ; nprnc(16)=nprnc(18) ; occ(16,1)=occ(18,1)
       kappa(17)= kappa(19) ; nprnc(17)=nprnc(19) ; occ(17,1)=occ(19,1)
       kappa(18) = k_h(1) ; nprnc(18) =  n_h(1)  ; occ(18,1)= o_h(1)
       kappa(19) = k_h(2) ; nprnc(19) =  n_h(2)  ; occ(19,1)= o_h(2)

       IF (jspins.EQ.2) THEN
          o_h(1) = occ(15,jspins) ; o_h(2) = occ(16,jspins)
          occ(15,jspins) = occ(17,jspins) 
          occ(16,jspins) = occ(18,jspins) 
          occ(17,jspins) = occ(19,jspins) 
          occ(18,jspins) = o_h(1) 
          occ(19,jspins) = o_h(2) 
       ENDIF
    ENDIF

  END SUBROUTINE setcor
      END MODULE m_setcor
