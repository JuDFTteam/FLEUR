!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_flipcdn
!     *******************************************************
!     this subroutine reads the charge density and flips the
!     magnetic moment within the m.t.sphere for each atom
!     according to the variable nflip. This variable is read in
!     the main program
!  TODO; (Test)           nflip = -1 : flip spin in sphere
!  TODO:           nflip = -2 : scale spin by bmu(n)
!             nflip = any: no spin flip
!                            r.pentcheva,kfa,Feb'96
!
!     Extension to multiple U per atom type by G.M. 2017
!
!     Removed integer nflip switch and added angles phi/theta
!     (and an additional spin scale switch)
!     which defines spin flip for each atom individually.
!     => Magnetisation axis can now be chosen independet
!     of spin quantization axis.
!     R. Hilgers, Okt. 2019
!     *******************************************************
CONTAINS

SUBROUTINE flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,phi,theta,optDen)
   USE m_rotdenmat
   USE m_constants
   USE m_cdn_io
   USE m_types

   IMPLICIT NONE

   TYPE(t_stars),INTENT(IN)    :: stars
   TYPE(t_vacuum),INTENT(IN)   :: vacuum
   TYPE(t_atoms),INTENT(IN)    :: atoms
   TYPE(t_sphhar),INTENT(IN)   :: sphhar
   TYPE(t_input),INTENT(IN)    :: input
   TYPE(t_sym),INTENT(IN)      :: sym
   TYPE(t_noco),INTENT(IN)     :: noco
   TYPE(t_oneD),INTENT(IN)     :: oneD
   TYPE(t_cell),INTENT(IN)     :: cell
   REAL, OPTIONAL, INTENT(IN)  :: phi(atoms%ntype)
   REAL, OPTIONAL, INTENT(IN)  :: theta(atoms%ntype)
   TYPE(t_potden), OPTIONAL,INTENT(INOUT) :: optDen

   ! Local type instance
   TYPE(t_potden)            :: den

   ! Local Scalars
   COMPLEX                   :: rhodummy, imPart12, realPart12
   REAL                      :: rhodumms,fermiEnergyTemp, realPart1, realPart2, imPart1,imPart2, rhodummyR, rotAnglePhi(atoms%ntype),rotAngleTheta(atoms%ntype)
   INTEGER                   :: i,nt,j,lh,na,mp,ispin,urec,itype,m,i_u,k
   INTEGER                   :: archiveType
   LOGICAL                   :: n_exist,l_qfix,l_error, l_flip(atoms%ntype), scaleSpin(atoms%ntype)
   ! Local Arrays
   CHARACTER(len=80), ALLOCATABLE :: clines(:)



!Flipcdn by optional given angle if lflip is false but routine is called.
   DO k=1, atoms%ntype
      IF(.NOT.input%lflip.AND.PRESENT(phi).AND.present(theta)) THEN
         rotAnglePhi(k)=phi(k)
         rotAngleTheta(k)=theta(k)
         scaleSpin(k)=.FALSE.
      ELSE IF (input%lflip) THEN
!Rotation triggerd by lflip.
         rotAnglePhi(k)=atoms%flipSpinPhi(k)
         rotAngleTheta(k)=atoms%flipSpinTheta(k)
         scaleSpin(k)=atoms%flipSpinScale(k)
      ELSE
         CALL judft_error("You shouldn't be here. There went something wrong.",calledby="flipcdn")
      END IF
   END DO



   DO itype=1, atoms%ntype
      l_flip(itype)=MERGE(.TRUE.,.FALSE.,(rotAnglePhi(itype).NE.0.0) .OR.(rotAngleTheta(itype).NE.0.0))
   END DO
   !rot_den_mat(alph,beta,rho11,rho22,rho21)
   IF(.NOT.PRESENT(optDen)) THEN
      CALL den%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
         IF(noco%l_noco) THEN
      archiveType = CDN_ARCHIVE_TYPE_NOCO_const
   ELSE
      archiveType = CDN_ARCHIVE_TYPE_CDN1_const
   END IF
      ! read the charge density
       CALL readDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,&
                    CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den)
   ELSE
      den=optDen
   END IF


   ! flip cdn for each atom with rotation angles given
   na = 1
   DO itype = 1, atoms%ntype
      IF (l_flip(itype).AND.(.NOT.scaleSpin(itype))) THEN
         ! spherical and non-spherical m.t. charge density
!$OMP parallel private(rhodummy,j,rhodummyR)
!$OMP do
	DO lh = 0,sphhar%nlh(sym%ntypsy(na))
            DO j = 1,atoms%jri(itype)
                IF (noco%l_mtNocoPot) THEN
                   rhodummy=CMPLX(den%mt(j,lh,itype,3),den%mt(j,lh,itype,4))
                   CALL rot_den_mat(rotAnglePhi(itype),rotAngleTheta(itype),den%mt(j,lh,itype,1),den%mt(j,lh,itype,2),rhodummy)
                   den%mt(j,lh,itype,3)=REAL(rhodummy)
		   den%mt(j,lh,itype,4)=AIMAG(rhodummy)
                ELSE
                   IF (rotAngleTheta(itype).EQ.(pimach()).AND.rotAnglePhi(itype).EQ.0) THEN
                      rhodummyR = den%mt(j,lh,itype,1)
                      den%mt(j,lh,itype,1) = den%mt(j,lh,itype,input%jspins)
                      den%mt(j,lh,itype,input%jspins) = rhodummyR
                   ELSE
                      !Since in non-noco case the den-matrices are only initialized with two diagonal components we cannot perform flips where off-diagonal elements arise in non-noco case => Only rotations by theta=Pi/2 are allowed.
                      CALL judft_error("l_mtNocoPot=F in combination with spin flips different from flipSpinTheta=Pi and flipSpinPhi=0 is currently not supported.",&
                         calledby="flipcdn")
                   END IF
                END IF
            END DO
         END DO
!$OMP end do
!$OMP end parallel
      ELSE IF (l_flip(itype).AND.scaleSpin(itype)) THEN
         IF((rotAngleTheta(itype).NE.(pimach()) .OR.rotAnglePhi(itype).NE.0.0)) CALL judft_error("Spinscaling in combination with flipSpin is currently only implemented using flipSpinTheta=Pi and flipSpinPhi=0.0.",calledby="flipcdn")
         DO lh = 0,sphhar%nlh(sym%ntypsy(na))
            DO j = 1,atoms%jri(itype)
               rhodummy = den%mt(j,lh,itype,1) + den%mt(j,lh,itype,input%jspins)
               rhodumms = den%mt(j,lh,itype,1) - den%mt(j,lh,itype,input%jspins)
               den%mt(j,lh,itype,1) = 0.5 * (rhodummy + atoms%bmu(itype)*rhodumms)
               den%mt(j,lh,itype,input%jspins) = 0.5 * (rhodummy - atoms%bmu(itype)*rhodumms )
            END DO
         END DO
      END IF
         na = na + atoms%neq(itype)
   END DO

   IF (input%l_onlyMtStDen) THEN
!!This Segment takes care that no interstitial magnetization is written in the the density. Meaning: Off diagonal elements of density matrix set to 0 and diagonal elements of density matrix are equal to their mean value.
      den%pw(:,2)=(den%pw(:,1)+den%pw(:,2))*0.5 !mean value
      den%pw(:,1)=den%pw(:,2)
      IF (noco%l_noco) THEN
         den%pw(:,3)=CMPLX(0.0,0.0)
      END IF
   END IF


   ! for LDA+U: flip density matrix
   IF (ANY(den%mmpMat(:,:,:,:).NE.0.0).AND.atoms%n_u+atoms%n_hia>0) THEN
      DO i_u = 1, atoms%n_u+atoms%n_hia
         itype = atoms%lda_u(i_u)%atomType
        IF (l_flip(itype).AND.(.NOT.scaleSpin(itype))) THEN
            DO m = -3,3
               DO mp = -3,3
                  IF (noco%l_mtNocoPot) THEN
! Since den%mmpMat is complex but rot_den_mat can only handle real values as diagonals of den_mat a splitting of den%mmpMat in real and imaginary part is performed. Rotations are performed seperatly and added up afterwards.
                    realPart1=REAL(den%mmpMat(m,mp,i_u,1))
                    realPart2=REAL(den%mmpMat(m,mp,i_u,2))
                    realPart12=CMPLX(REAL(den%mmpMat(m,mp,i_u,3)),0)
                    imPart1=AIMAG(den%mmpMat(m,mp,i_u,1))
                    imPart2=AIMAG(den%mmpMat(m,mp,i_u,2))
                    imPart12=CMPLX(0,AIMAG(den%mmpMat(m,mp,i_u,3)))
                    CALL rot_den_mat(rotAnglePhi(itype),rotAngleTheta(itype),realPart1,realPart2,&
                      realPart12)
                    CALL rot_den_mat(rotAnglePhi(itype),rotAngleTheta(itype),imPart1,imPart2,&
                      imPart12)
                    den%mmpMat(m,mp,i_u,1)=CMPLX(realPart1,imPart1)
                    den%mmpMat(m,mp,i_u,2)=CMPLX(realPart2,imPart2)
                    den%mmpMat(m,mp,i_u,3)=realPart12+imPart12
                 ELSE
                    IF (rotAngleTheta(itype).EQ.(pimach()).AND.rotAnglePhi(itype).EQ.0) THEN
                       rhodummyR = den%mmpMat(m,mp,i_u,1)
                       den%mmpMat(m,mp,i_u,1) = den%mmpMat(m,mp,i_u,input%jspins)
                       den%mmpMat(m,mp,i_u,input%jspins) = rhodummyR
                    ELSE
                       !Since in non-noco case the den-matrices are only initialized with two diagonal components we cannot perform flips where off-diagonal elements arise in non-noco case => Only rotations by Pi degrees are allowed.
                       CALL judft_error("l_mtNocoPot=F in combination with spin flips different from flipSpinTheta=Pi, flipSpinPhi=0 is currently not supported.",&
                         calledby="flipcdn")
                    END IF
                 END IF

END DO
            END DO
          ELSE IF (l_flip(itype).AND.(scaleSpin(itype))) THEN
            DO m = -3,3
               DO mp = -3,3
                  IF((rotAngleTheta(itype).NE.pimach() .OR.rotAnglePhi(itype).NE.0.0)) CALL judft_error("Spinscaling in combination with flipSpin is currently only implemented using flipSpinTheta=Pi and flipSpinPhi=0.0",calledby="flipcdn")
                  rhodummy = den%mmpMat(m,mp,i_u,1) + den%mmpMat(m,mp,i_u,input%jspins)
                  rhodumms = den%mmpMat(m,mp,i_u,1) - den%mmpMat(m,mp,i_u,input%jspins)
                  den%mmpMat(m,mp,i_u,1) = 0.5 * (rhodummy + atoms%bmu(itype) * rhodumms)
                  den%mmpMat(m,mp,i_u,input%jspins) = 0.5 * (rhodummy - atoms%bmu(itype) * rhodumms)
               END DO
            END DO
         END IF
      END DO
   END IF

   ! write the spin-polarized density
   CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                     0,-1.0,0.0,.FALSE.,den)

   ! read enpara and  flip lines
   INQUIRE(file='enpara',exist=n_exist)
   IF (n_exist) THEN
      OPEN(40,file ='enpara',status='old',form='formatted')

      j = 2
      DO itype = 1, atoms%ntype
         j = j + 1
         IF (atoms%nlo(itype)>0) j = j + 2
      END DO
      IF (input%film) j = j + 1
      ALLOCATE (clines(2*j))
      DO i = 1, 2*j
         READ (40,'(a)') clines(i)
      END DO

      REWIND 40
      i = 0
      DO ispin = 1,input%jspins
         i = i + 2
         WRITE (40,'(a)') TRIM(clines(i-1))
         WRITE (40,'(a)') TRIM(clines(i))
         DO itype = 1, atoms%ntype
            i = i + 1
            m = i
            IF (l_flip(itype)) m = MOD(i+j,2*j)
            IF (m==0) m = 2*j
            WRITE (40,'(a)') TRIM(clines(m))
            IF (atoms%nlo(itype)>0) THEN
               WRITE (40,'(a)') TRIM(clines(m+1))
               WRITE (40,'(a)') TRIM(clines(m+2))
               i = i + 2
            END IF
         END DO
         IF (input%film) THEN
            i = i + 1
            WRITE (40,'(a)') TRIM(clines(i))
         END IF
      END DO
      DEALLOCATE (clines)
      CLOSE(40)
   END IF


IF(PRESENT(optDen)) optDen=den

END SUBROUTINE flipcdn

END MODULE m_flipcdn
