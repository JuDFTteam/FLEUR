!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_coredr
   implicit none
CONTAINS
  SUBROUTINE coredr(input,atoms,iType,seig, rho,sphhar, vrs, qints,rhc,l_useOtherCoreSolver)
    !     *******************************************************
    !     *****   set up the core densities for compounds   *****
    !     *****   for relativistic core                     *****
    !     *******************************************************

    USE m_etabinit
    USE m_spratm
    USE m_ccdnup
    USE m_cdn_io
    USE m_types_input
    USE m_types_atoms
    USE m_types_sphhar
    
    IMPLICIT NONE
    
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    !
    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: iType
    REAL,    INTENT(OUT) :: seig
    !     ..
    !     .. Array Arguments ..
    REAL   , INTENT (IN)    :: vrs(atoms%jmtd,atoms%ntype,input%jspins)
    REAL,    INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
    REAL,    INTENT (OUT)   :: rhc(atoms%msh,atoms%ntype,input%jspins),qints(atoms%ntype,input%jspins)
    LOGICAL, INTENT (INOUT) :: l_useOtherCoreSolver
    !     ..
    !     .. Local Scalars ..
    REAL dxx,rnot,sume,t2,t2b,z,t1,rr,d,v1,v2
    INTEGER i,j,jspin,k,ncmsh
    !     ..
    !     .. Local Arrays ..
    REAL br(atoms%jmtd,atoms%ntype),brd(atoms%msh),etab(100,atoms%ntype),&
         rhcs(atoms%jmtd,atoms%ntype,input%jspins),rhochr(atoms%msh),rhospn(atoms%msh),&
         tecs(atoms%ntype,input%jspins),vr(atoms%jmtd,atoms%ntype),vrd(atoms%msh)
    INTEGER nkmust(atoms%ntype),ntab(100,atoms%ntype),ltab(100,atoms%ntype)

    !     ..
    ntab(:,:) = -1 ; ltab(:,:) = -1 ; etab(:,:) = 0.0
    !
    ! setup potential and field
    !
    
    seig = 0.0
    
    IF (input%jspins.EQ.1) THEN
       DO j = 1,atoms%jmtd
          vr(j,iType) = vrs(j,iType,1)
          br(j,iType) = 0.0
       END DO
    ELSE
       DO j = 1,atoms%jmtd
          vr(j,iType) = (vrs(j,iType,1)+vrs(j,iType,input%jspins))/2.
          br(j,iType) = (vrs(j,iType,input%jspins)-vrs(j,iType,1))/2.
       END DO
       IF(MAXVAL(ABS(br(1:atoms%jmtd,iType))).LT.1.0e-8) THEN
          l_useOtherCoreSolver = .TRUE. ! Use the other solver in case of a nonmagnetic atom in a magnetic calculation.
          RETURN
       END IF
    END IF
    !
    ! setup eigenvalues
    CALL etabinit(atoms,input, vr, etab,ntab,ltab,nkmust)
    !
    ncmsh = atoms%msh
    ! ---> set up densities

    DO j = 1,atoms%jri(iType)
       vrd(j) = vr(j,iType)
       brd(j) = br(j,iType)
    END DO

    IF (input%l_core_confpot) THEN
       !--->    linear extension of the potential with slope t1 / a.u.
       rr = atoms%rmt(iType)
       d = EXP(atoms%dx(iType))
       t1=0.125
       !         t2  = vrd(jri(iType))/rr - rr*t1
       !         t2b = brd(jri(iType))/rr - rr*t1
       t2  = vrs(atoms%jri(iType),iType,1)     /rr - rr*t1
       t2b = vrs(atoms%jri(iType),iType,input%jspins)/rr - rr*t1
    ELSE
       t2 = vrd(atoms%jri(iType))/ (atoms%jri(iType)-ncmsh)
       t2b = brd(atoms%jri(iType))/ (atoms%jri(iType)-ncmsh)
    ENDIF
    IF (atoms%jri(iType).LT.ncmsh) THEN
       DO i = atoms%jri(iType) + 1,ncmsh
          IF (input%l_core_confpot) THEN
             rr = d*rr
             v1 = rr*( t2  + rr*t1 )
             v2 = rr*( t2b + rr*t1 )
             vrd(i) = 0.5*(v2 + v1)
             brd(i) = 0.5*(v2 - v1)
          ELSE
             vrd(i) = vrd(atoms%jri(iType)) + t2* (i-atoms%jri(iType))
             brd(i) = brd(atoms%jri(iType)) + t2b* (i-atoms%jri(iType))
          ENDIF
       END DO
    END IF

    !
    rnot = atoms%rmsh(1,iType)
    z = atoms%zatom(iType)
    dxx = atoms%dx(iType)

    CALL spratm(atoms%msh,vrd,brd,z,rnot,dxx,ncmsh,etab(1,iType),ntab(1,iType),ltab(1,iType), sume,rhochr,rhospn)

    seig = seig + atoms%neq(iType)*sume
    !
    !     rho_up=2(ir) = (rhochr(ir)  + rhospn(ir))*0.5
    !     rho_dw=1(ir) = (rhochr(ir)  - rhospn(ir))*0.5
    !
    IF (input%jspins.EQ.2) THEN
       DO j = 1,atoms%jri(iType)
          rhcs(j,iType,input%jspins) = (rhochr(j)+rhospn(j))*0.5
          rhcs(j,iType,1) = (rhochr(j)-rhospn(j))*0.5
       END DO
    ELSE
       DO j = 1,atoms%jri(iType)
          rhcs(j,iType,1) = rhochr(j)
       END DO
    END IF
    IF (input%jspins.EQ.2) THEN
       DO j = 1,atoms%msh
          rhc(j,iType,input%jspins) = (rhochr(j)+rhospn(j))*0.5
          rhc(j,iType,1) = (rhochr(j)-rhospn(j))*0.5
       ENDDO
    ELSE
       DO j = 1,atoms%msh
          rhc(j,iType,1) = rhochr(j)
       END DO
    END IF
    !
    !---->update spherical charge density rho with the core density.
    CALL ccdnup(atoms,sphhar,input,iType, rho, sume,vrs,rhochr,rhospn, tecs,qints)
    !
    !----> store core charge densities
    CALL writeCoreDensity(input,atoms,rhcs,tecs,qints)
  END SUBROUTINE coredr
END MODULE m_coredr
