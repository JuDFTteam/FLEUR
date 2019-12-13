MODULE m_analysistests
   USE m_types

   INTEGER, PARAMETER :: numtests=3
   INTEGER, PARAMETER :: numders=5
   ! 1 = grdchlh with ndvgrd=3
   ! 5 = SPEX Derivative
   INTEGER, PARAMETER :: numints=5

CONTAINS

   SUBROUTINE buildfunctions(atoms,radii,funcsr,trueder,trueint)

      TYPE(t_atoms), INTENT(IN) :: atoms
      REAL, INTENT(OUT) :: radii(atoms%jmtd,atoms%ntype), funcsr(atoms%jmtd,atoms%ntype,numtests), trueder(atoms%jmtd,atoms%ntype,numtests), trueint(atoms%jmtd,atoms%ntype,numtests)
     
      REAL :: Rmt(atoms%jmtd,atoms%ntype)
      INTEGER :: n
      
      radii=atoms%rmsh
      DO n=1,atoms%ntype
         Rmt(:,n)=atoms%rmt(n)
      END DO

      funcsr(:,:,1)=radii(:,:)**2
      trueder(:,:,1)=2*radii(:,:)
      trueint(:,:,1)=(1./3.)*radii(:,:)**3

      funcsr(:,:,2)=radii(:,:)**3
      trueder(:,:,2)=3*radii(:,:)**2
      trueint(:,:,2)=(1./4.)*radii(:,:)**4

      funcsr(:,:,3)=radii(:,:)**2-(radii(:,:)**3)/Rmt(:,:)
      trueder(:,:,3)=2*radii(:,:)-3*(radii(:,:)**2)/Rmt(:,:)
      trueint(:,:,3)=(1./3.)*radii(:,:)**3-(1./4.)*(radii(:,:)**4)/Rmt(:,:)

   END SUBROUTINE buildfunctions
   
   SUBROUTINE derivtest(atoms,radii,funcsr,testders)
      USE m_grdchlh
      USE m_gradYlm

      TYPE(t_atoms), INTENT(IN) :: atoms
      REAL, INTENT(IN) :: radii(atoms%jmtd,atoms%ntype), funcsr(atoms%jmtd,atoms%ntype,numtests)
      REAL, INTENT(OUT) :: testders(atoms%jmtd,atoms%ntype,numtests,numders)

      INTEGER :: n
      REAL :: junkder(atoms%jmtd)

      DO n=1,atoms%ntype
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,1),3, testders(:,n,1,1),junkder)
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,1),4, testders(:,n,1,2),junkder)
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,1),5, testders(:,n,1,3),junkder)
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,1),6, testders(:,n,1,4),junkder)
         CALL Derivative(funcsr(:,n,1), n, atoms, testders(:,n,1,5))
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,2),3, testders(:,n,2,1),junkder)
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,2),4, testders(:,n,2,2),junkder)
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,2),5, testders(:,n,2,3),junkder)
         CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),radii(:,n),funcsr(:,n,2),6, testders(:,n,2,4),junkder)
         CALL Derivative(funcsr(:,n,2), n, atoms, testders(:,n,2,5))
         CALL Derivative(funcsr(:,n,3), n, atoms, testders(:,n,3,5))
      END DO
   
   END SUBROUTINE derivtest

   SUBROUTINE integtest(atoms,radii,funcsr,testints)
      USE m_intgr

      TYPE(t_atoms), INTENT(IN) :: atoms
      REAL, INTENT(IN) :: radii(atoms%jmtd,atoms%ntype), funcsr(atoms%jmtd,atoms%ntype,numtests)
      REAL, INTENT(OUT) :: testints(atoms%jmtd,atoms%ntype,numtests,numints)

      INTEGER :: n

      DO n=1,atoms%ntype      
         CALL intgr2(funcsr(:atoms%jri(n),n,1),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,1,1))
         CALL intgr2(funcsr(:atoms%jri(n),n,2),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,2,1))
         CALL intgr2(funcsr(:atoms%jri(n),n,3),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,3,1))
         CALL intgrt(funcsr(:atoms%jri(n),n,1),radii(:,n),atoms%jri(n),testints(:atoms%jri(n),n,1,2))
         CALL intgrt(funcsr(:atoms%jri(n),n,2),radii(:,n),atoms%jri(n),testints(:atoms%jri(n),n,2,2))
         CALL intgrt(funcsr(:atoms%jri(n),n,3),radii(:,n),atoms%jri(n),testints(:atoms%jri(n),n,3,2))
         CALL intgrtlog(funcsr(:atoms%jri(n),n,1),radii(:,n),atoms%jri(n),testints(:atoms%jri(n),n,1,3))
         CALL intgrtlog(funcsr(:atoms%jri(n),n,2),radii(:,n),atoms%jri(n),testints(:atoms%jri(n),n,2,3))
         CALL intgrtlog(funcsr(:atoms%jri(n),n,3),radii(:,n),atoms%jri(n),testints(:atoms%jri(n),n,3,3))
         CALL intgr4(funcsr(:atoms%jri(n),n,1),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,1,4))
         CALL intgr4(funcsr(:atoms%jri(n),n,2),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,2,4))
         CALL intgr4(funcsr(:atoms%jri(n),n,3),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,3,4))
         CALL intgr5(funcsr(:atoms%jri(n),n,1),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,1,5))
         CALL intgr5(funcsr(:atoms%jri(n),n,2),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,2,5))
         CALL intgr5(funcsr(:atoms%jri(n),n,3),radii(:,n),atoms%dx(n),atoms%jri(n),testints(:atoms%jri(n),n,3,5))
      END DO
   
   END SUBROUTINE integtest

END MODULE m_analysistests
