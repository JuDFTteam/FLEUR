MODULE m_coredr
CONTAINS
  SUBROUTINE coredr(input,atoms,seig, rho,DIMENSION,sphhar, vrs, qints,rhc)
    !     *******************************************************
    !     *****   set up the core densities for compounds   *****
    !     *****   for relativistic core                     *****
    !     *******************************************************

    USE m_etabinit
    USE m_spratm
    USE m_ccdnup

    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    !
    !     .. Scalar Arguments ..
    REAL seig
    !     ..
    !     .. Array Arguments ..
    REAL   , INTENT (IN) :: vrs(atoms%jmtd,atoms%ntype,DIMENSION%jspd)
    REAL,    INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd)
    REAL,    INTENT (OUT) :: rhc(DIMENSION%msh,atoms%ntype,DIMENSION%jspd),qints(atoms%ntype,DIMENSION%jspd)
    !     ..
    !     .. Local Scalars ..
    REAL dxx,rnot,sume,t2,t2b,z,t1,rr,d,v1,v2
    INTEGER i,j,jatom,jspin,k,n,ncmsh
    LOGICAL exetab
    !     ..
    !     .. Local Arrays ..
    REAL br(atoms%jmtd,atoms%ntype),brd(DIMENSION%msh),etab(100,atoms%ntype),&
         rhcs(atoms%jmtd,atoms%ntype,DIMENSION%jspd),rhochr(DIMENSION%msh),rhospn(DIMENSION%msh),&
         tecs(atoms%ntype,DIMENSION%jspd),vr(atoms%jmtd,atoms%ntype),vrd(DIMENSION%msh)
    INTEGER nkmust(atoms%ntype),ntab(100,atoms%ntype),ltab(100,atoms%ntype)

    !     ..
    ntab(:,:) = -1 ; ltab(:,:) = -1 ; etab(:,:) = 0.0
    !
    ! setup potential and field
    !
    IF (input%jspins.EQ.1) THEN
       DO n = 1,atoms%ntype
          DO j = 1,atoms%jmtd
             vr(j,n) = vrs(j,n,1)
             br(j,n) = 0.0
          END DO
       END DO
    ELSE
       DO n = 1,atoms%ntype
          DO j = 1,atoms%jmtd
             vr(j,n) = (vrs(j,n,1)+vrs(j,n,input%jspins))/2.
             br(j,n) = (vrs(j,n,input%jspins)-vrs(j,n,1))/2.
          END DO
       END DO
    END IF
    !
    ! setup eigenvalues
    exetab = .FALSE.
    INQUIRE (file='core.dat',exist=exetab)
    IF (exetab) THEN
       OPEN (58,file='core.dat',form='formatted',status='old')
       REWIND 58
       DO n = 1,atoms%ntype
          READ (58,FMT=*) nkmust(n)
          DO k = 1,nkmust(n)
             READ (58,FMT='(f12.6,2i3)') etab(k,n),ntab(k,n),&
                  &                                               ltab(k,n)

          END DO
       END DO
    ELSE
       OPEN (58,file='core.dat',form='formatted',status='new')
       CALL etabinit(atoms,DIMENSION,input, vr, etab,ntab,ltab,nkmust)
    END IF
    !
    ncmsh = DIMENSION%msh
    seig = 0.
    ! ---> set up densities
    DO jatom = 1,atoms%ntype
       !
       DO j = 1,atoms%jri(jatom)
          vrd(j) = vr(j,jatom)
          brd(j) = br(j,jatom)
       END DO

       IF (input%l_core_confpot) THEN
          !--->    linear extension of the potential with slope t1 / a.u.
          rr = atoms%rmt(jatom)
          d = EXP(atoms%dx(jatom))
          t1=0.125
          !         t2  = vrd(jri(jatom))/rr - rr*t1
          !         t2b = brd(jri(jatom))/rr - rr*t1
          t2  = vrs(atoms%jri(jatom),jatom,1)     /rr - rr*t1
          t2b = vrs(atoms%jri(jatom),jatom,input%jspins)/rr - rr*t1
       ELSE
          t2 = vrd(atoms%jri(jatom))/ (atoms%jri(jatom)-ncmsh)
          t2b = brd(atoms%jri(jatom))/ (atoms%jri(jatom)-ncmsh)
       ENDIF
       IF (atoms%jri(jatom).LT.ncmsh) THEN
          DO i = atoms%jri(jatom) + 1,ncmsh
             IF (input%l_core_confpot) THEN
                rr = d*rr
                v1 = rr*( t2  + rr*t1 )
                v2 = rr*( t2b + rr*t1 )
                vrd(i) = 0.5*(v2 + v1)
                brd(i) = 0.5*(v2 - v1)
             ELSE
                vrd(i) = vrd(atoms%jri(jatom)) + t2* (i-atoms%jri(jatom))
                brd(i) = brd(atoms%jri(jatom)) + t2b* (i-atoms%jri(jatom))
             ENDIF
          END DO
       END IF

       !        rr = rmsh(1,jatom)
       !        do i =1, ncmsh
       !          rr = d*rr
       !         write(*,'(3f20.10)') rr,vrd(i),brd(i)
       !        enddo

       !
       rnot = atoms%rmsh(1,jatom)
       z = atoms%zatom(jatom)
       dxx = atoms%dx(jatom)

       CALL spratm(DIMENSION%msh,vrd,brd,z,rnot,dxx,ncmsh,&
            etab(1,jatom),ntab(1,jatom),ltab(1,jatom), sume,rhochr,rhospn)

       seig = seig + atoms%neq(jatom)*sume
       !
       !     rho_up=2(ir) = (rhochr(ir)  + rhospn(ir))*0.5
       !     rho_dw=1(ir) = (rhochr(ir)  - rhospn(ir))*0.5
       !
       IF (input%jspins.EQ.2) THEN
          DO j = 1,atoms%jri(jatom)
             rhcs(j,jatom,input%jspins) = (rhochr(j)+rhospn(j))*0.5
             rhcs(j,jatom,1) = (rhochr(j)-rhospn(j))*0.5
          END DO
       ELSE
          DO j = 1,atoms%jri(jatom)
             rhcs(j,jatom,1) = rhochr(j)
          END DO
       END IF
       IF (input%jspins.EQ.2) THEN
          DO j = 1,DIMENSION%msh
             rhc(j,jatom,input%jspins) = (rhochr(j)+rhospn(j))*0.5
             rhc(j,jatom,1) = (rhochr(j)-rhospn(j))*0.5
          ENDDO
       ELSE
          DO j = 1,DIMENSION%msh
             rhc(j,jatom,1) = rhochr(j)
          END DO
       END IF
       !
       ! store atomic eigenvalues to file.58
       IF (jatom.EQ.1) REWIND 58
       WRITE (58,FMT=*) nkmust(jatom)
       DO k = 1,nkmust(jatom)
          WRITE (58,FMT='(f12.6,2i3)') etab(k,jatom),ntab(k,jatom), ltab(k,jatom)
       END DO
       !---->update spherical charge density rho with the core density.
       CALL ccdnup(atoms,sphhar,input,jatom, rho, sume,vrs,rhochr,rhospn, tecs,qints)

    END DO ! loop over atoms (jatom)
    !
    !---->store core charge densities to file.17
    OPEN (17,file='cdnc',form='unformatted',status='unknown')
    REWIND 17
    DO jspin = 1,input%jspins
       DO jatom = 1,atoms%ntype
          WRITE (17) (rhcs(j,jatom,jspin),j=1,atoms%jri(jatom))
          WRITE (17) tecs(jatom,jspin)
       END DO
       WRITE (17) (qints(jatom,jspin),jatom=1,atoms%ntype)
    END DO
    CLOSE (17)
    !
    RETURN
  END SUBROUTINE coredr
END MODULE m_coredr
