MODULE m_types_profile

   TYPE :: t_profile
      REAL :: kmax
      REAL :: rmtFactor
      REAL :: lmaxFactor
      REAL :: fermiSmearing

      CHARACTER(LEN=20) :: profileName
      CHARACTER(LEN=20) :: addLOSetup

      CONTAINS

      PROCEDURE :: init => initProfile
      PROCEDURE :: load => loadProfile
   END TYPE t_profile

   PUBLIC :: t_profile

   CONTAINS

   SUBROUTINE initProfile(this)

      IMPLICIT NONE

      CLASS(t_profile), INTENT(INOUT) :: this

      this%profileName = "default"
      this%kmax = 4.5
      this%rmtFactor = 1.0
      this%lmaxFactor = 1.0
      this%fermiSmearing = 0.001
      this%addLOSetup = ""

   END SUBROUTINE initProfile

   SUBROUTINE loadProfile(this,profileName)

      IMPLICIT NONE

      CLASS(t_profile), INTENT(INOUT) :: this

      CHARACTER(LEN=*), INTENT(IN) :: profileName


      INTEGER :: io_stat
      REAL    :: kmax, rmtFactor, lmaxFactor, fermiSmearing
      LOGICAL :: l_exist, l_found

      CHARACTER(LEN=20) :: name
      CHARACTER(LEN=20) :: addLOSetup
      CHARACTER(LEN=20) :: filename
      CHARACTER(len=8)  :: str

      NAMELIST /profile/ name,kmax,rmtFactor,lmaxFactor,addLOSetup,fermiSmearing

      filename = "profile.config"

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
            rmtFactor = 1.0
            lmaxFactor = 1.0
            addLOSetup = ""
            fermiSmearing = 0.001
            READ(558,profile,iostat=io_stat)
            IF (io_stat.EQ.0) THEN
               IF(TRIM(ADJUSTL(profileName)).EQ.TRIM(ADJUSTL(name))) THEN
                  this%profileName = name
                  this%kmax = kmax
                  this%rmtFactor = rmtFactor
                  this%lmaxFactor = lmaxFactor
                  this%addLOSetup = TRIM(ADJUSTL(addLOSetup))
                  this%fermiSmearing = fermiSmearing
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
