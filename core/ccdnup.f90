MODULE m_ccdnup
  !     *******************************************************
  !     *****   set up the core densities for compounds.  *****
  !     *****   in accordanse to d.d.koelling's cored     *****
  !     *******************************************************
CONTAINS
  SUBROUTINE ccdnup(&
       &                  atoms,sphhar,input,jatom,&
       &                  rho,&
       &                  sume,vrs,rhochr,rhospn,&
       &                  tecs,qints)

    USE m_constants,ONLY:sfp_const
    USE m_intgr, ONLY : intgr3
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jatom 
    REAL,    INTENT (IN) :: sume
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rhochr(:),rhospn(:)!(dimension%msh)
    REAL,    INTENT (IN) :: vrs(:,:,:)!(atoms%jmtd,atoms%ntype,input%jspins)
    REAL,    INTENT (OUT) :: tecs(:,:)!(atoms%ntype,input%jspins)
    REAL,    INTENT (OUT) :: qints(:,:)!(atoms%ntype,input%jspins)
    REAL,    INTENT (INOUT) :: rho(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
    !     ..
    !     .. Local Scalars ..
    REAL d,dxx,q,rad,rhs
    INTEGER i,j,jspin,nm,nm1
    !     ..
    !     .. Local Arrays ..
    REAL rhoc(SIZE(rhochr)),rhoss(SIZE(rhochr))
    !     ..
    nm = atoms%jri(jatom)
    !     ---->update spherical charge density rho with the core density.
    !     ---->for spin-polarized (jspins=2), take only half the density
    DO jspin = 1,input%jspins
       IF (input%jspins.EQ.2 .AND. jspin.EQ.1) THEN
          DO j = 1,SIZE(rhochr)
             rhoss(j) = rhochr(j) - rhospn(j)
          END DO
       ELSE IF (input%jspins.EQ.2 .AND. jspin.EQ.2) THEN
          DO j = 1,SIZE(rhochr)
             rhoss(j) = rhochr(j) + rhospn(j)
          END DO
          ! jspins=1
       ELSE
          DO j = 1,SIZE(rhochr)
             rhoss(j) = rhochr(j)
          END DO
          !
       END IF
       !
       DO  j = 1,nm
          rhoc(j) = rhoss(j)/input%jspins
          rho(j,0,jatom,jspin) = rho(j,0,jatom,jspin) + rhoc(j)/sfp_const
       ENDDO
       DO  i = 1,nm
          rhoc(i) = rhoc(i)*vrs(i,jatom,jspin)/atoms%rmsh(i,jatom)
       ENDDO
       CALL intgr3(rhoc,atoms%rmsh(1,jatom),atoms%dx(jatom),nm,rhs)
       tecs(jatom,jspin) = sume/input%jspins - rhs
       WRITE (6,FMT=8010) jatom,jspin,tecs(jatom,jspin),sume/input%jspins
  
       !     ---> simpson integration
       dxx = atoms%dx(jatom)
       d = EXP(atoms%dx(jatom))
       rad = atoms%rmt(jatom)
       q = rad*rhoss(nm)/2.
       DO  nm1 = nm + 1,SIZE(rhochr) - 1,2
          rad = d*rad
          q = q + 2*rad*rhoss(nm1)
          rad = d*rad
          q = q + rad*rhoss(nm1+1)
       ENDDO
       q = 2*q*dxx/3
       WRITE (6,FMT=8000) q/input%jspins
       qints(jatom,jspin) = q*atoms%neq(jatom)

    END DO ! end-do-loop input%jspins

8000 FORMAT (f20.8,' electrons lost from core.')
8010 FORMAT (10x,'atom type',i3,'  (spin',i2,') ',/,10x,&
         &       'kinetic energy=',e20.12,5x,'sum of the eigenvalues=',&
         &       e20.12)

  END SUBROUTINE ccdnup
END MODULE m_ccdnup
