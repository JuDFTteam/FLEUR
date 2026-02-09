!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_profile
   implicit none
   private
   INTERFACE
      FUNCTION dropDefaultEConfig() BIND(C, name="dropDefaultEconfig")
      USE iso_c_binding
      INTEGER(c_int) dropDefaultEConfig
      END FUNCTION dropDefaultEConfig

      FUNCTION dropDefault2EConfig() BIND(C, name="dropDefault2EConfig")
      USE iso_c_binding
      INTEGER(c_int) dropDefault2EConfig
      END FUNCTION dropDefault2EConfig

      FUNCTION dropOxidesValEConfig() BIND(C, name="dropOxidesValidationEConfig")
      USE iso_c_binding
      INTEGER(c_int) dropOxidesValEConfig
      END FUNCTION dropOxidesValEConfig

      FUNCTION dropProfiles() BIND(C, name="dropProfiles")
      USE iso_c_binding
      INTEGER(c_int) dropProfiles
      END FUNCTION dropProfiles
   END INTERFACE



   TYPE :: t_profile
      REAL :: kmax ! This is K_max
      REAL :: kGmaxFactor ! G_max = G_maxXC = K_max * kGmaxFactor
      REAL :: rmtFactor ! This is a postprocessing factor to reduce the MT radii after their initial calculation
      REAL :: lmaxFactor ! lmax = Kmax * R_MT * lmaxfactor
      REAL :: fermiSmearing ! The Fermi smearing energy
      REAL :: kPDen ! The k-Point density

      CHARACTER(LEN=20) :: profileName
      CHARACTER(LEN=50) :: addLOSetup
      CHARACTER(LEN=20) :: atomSetup

      CONTAINS

      PROCEDURE :: init => initProfile
      PROCEDURE :: load => loadProfile
   END TYPE t_profile

   PUBLIC :: t_profile

   CONTAINS

   SUBROUTINE initProfile(this,profileName)

      IMPLICIT NONE

      CLASS(t_profile), INTENT(INOUT) :: this
      character(len=*), intent(in) :: profileName
      logical :: l_exist
      integer :: idum


     
      !Load profile from file if available
      INQUIRE(file='profile.config',exist=l_exist)
      IF (.NOT.l_exist) idum=dropProfiles()

      CALL this%load(profileName)

      !Drop atomic profiles if not existing
      IF(this%atomSetup.EQ."oxides_validation") THEN
         INQUIRE(file='oxides_validation.econfig',exist=l_exist)
         IF (.NOT.l_exist) idum=dropOxidesValEconfig()
      ELSE IF (this%atomSetup.EQ."default2") THEN
         INQUIRE(file='default2.econfig',exist=l_exist)
         IF (.NOT.l_exist) idum=dropDefault2Econfig()
      ELSE
         INQUIRE(file='default.econfig',exist=l_exist)
         IF (.NOT.l_exist) idum=dropDefaultEconfig()
      END IF

   END SUBROUTINE initProfile

   SUBROUTINE loadProfile(this,profileName_IN)

      IMPLICIT NONE

      CLASS(t_profile), INTENT(INOUT) :: this

      CHARACTER(LEN=*), INTENT(IN) :: profileName_IN


      INTEGER :: io_stat
      REAL    :: kmax, kGmaxFactor, rmtFactor, lmaxFactor, fermiSmearing, kPDen
      LOGICAL :: l_exist, l_found

      CHARACTER(LEN=20) :: name
      CHARACTER(LEN=50) :: addLOSetup
      CHARACTER(LEN=20) :: filename
      CHARACTER(len=8)  :: str
      CHARACTER(LEN=20) :: atomSetup
      CHARACTER(LEN=20) :: profileName

      NAMELIST /profile/ name,kmax,rmtFactor,lmaxFactor,addLOSetup,fermiSmearing,atomSetup,kGmaxFactor,kPDen

      filename = "profile.config"
      profileName = TRIM(ADJUSTL(profileName_IN))
      if (profileName.EQ."") profileName="default"

      INQUIRE(file=TRIM(filename),exist=l_exist)

      IF (.NOT.l_exist) THEN
         WRITE(*,*) 'No ', TRIM(filename), ' file found. Using default profile.'
         RETURN
      END IF

      l_found = .FALSE.

      OPEN(558,file=filename)
      DO
         READ(558,*,err=100,END=100) str
         IF (str.EQ."&profile") THEN
            BACKSPACE(558)
            name = "unknown"
            kmax = 4.5
            kGmaxFactor = 3.0
            rmtFactor = 1.0
            lmaxFactor = 1.0
            addLOSetup = ""
            fermiSmearing = 0.001
            kPDen = -1.0
            atomSetup = ""
            READ(558,profile,iostat=io_stat)
            IF (io_stat.EQ.0) THEN
               IF(profileName.EQ.TRIM(ADJUSTL(name))) THEN
                  this%profileName = name
                  this%kmax = kmax
                  this%kGmaxFactor = kGmaxFactor
                  this%rmtFactor = rmtFactor
                  this%lmaxFactor = lmaxFactor
                  this%addLOSetup = TRIM(ADJUSTL(addLOSetup))
                  this%fermiSmearing = fermiSmearing
                  this%kPDen = kPDen
                  this%atomSetup = TRIM(ADJUSTL(atomSetup))
                  l_found = .TRUE.
               END IF
            ELSE
               WRITE(*,*) 'Error in loading profile!'
               WRITE(*,*) 'io_stat: ', io_stat
            END IF
         END IF
      END DO
100   CLOSE(558)
      
      IF(.NOT.l_found) THEN
         WRITE(*,*) 'Could not find profile ', TRIM(ADJUSTL(profileName)), '. Using default profile.'
      ELSE
         WRITE(*,*) 'Using profile "', TRIM(ADJUSTL(profileName)), '".'
      END IF

   END SUBROUTINE loadProfile




END MODULE m_types_profile
