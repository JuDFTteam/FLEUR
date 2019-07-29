MODULE m_types_hub1ham

   USE m_constants

   IMPLICIT NONE

   PRIVATE

   TYPE t_hub1ham
      !Contains the information for the hubbard 1 solver
      !information for the atomic hamiltonian
      INTEGER                    :: iter
      LOGICAL                    :: l_runthisiter   !switch which determines wether Hubbard 1 will be run in the current iteration
      REAL, ALLOCATABLE          :: init_occ(:)     !initial occupation
      REAL                       :: beta            !inverse Temperature
      INTEGER                    :: n_exc = 2       !number of excitations considered 
 
      
      !Switches for args that were explicitly given and should not be calculated from DFT
      LOGICAL,ALLOCATABLE :: l_soc_given(:) 
      LOGICAL,ALLOCATABLE :: l_ccf_given(:) 

      !Additional arguments to be passed on to hloc.cfg (at the moment only real)
      INTEGER,             ALLOCATABLE :: n_addArgs(:) 
      CHARACTER(len=50),   ALLOCATABLE :: arg_keys(:,:)
      REAL,                ALLOCATABLE :: arg_vals(:,:)

      !Exchange splitting
      INTEGER, ALLOCATABLE       :: n_exc_given(:)
      REAL,    ALLOCATABLE       :: exc(:,:)        !exchange splitting
      INTEGER, ALLOCATABLE       :: exc_l(:,:)      !l quantum number from which the intraorbital exchange 
      REAL,    ALLOCATABLE       :: mag_mom(:,:)    !magnetic moment (for exchange splitting)
      REAL, ALLOCATABLE          :: init_mom(:,:)     !initial magnetic moment

      !These are the properties calculated from DFT
      REAL, ALLOCATABLE          :: xi(:)           !Spin-orbit coupling parameter
      REAL, ALLOCATABLE          :: ccf(:)          !crystal field factor
      REAL, ALLOCATABLE          :: ccfmat(:,:,:)   !crystal field splitting matrix

      CONTAINS 

      PROCEDURE, PASS :: init => hub1_init

   END TYPE t_hub1ham

   PUBLIC t_hub1ham

   CONTAINS 

   SUBROUTINE hub1_init(this,n_maxhub,n_maxargs)

      IMPLICIT NONE

      CLASS(t_hub1ham),    INTENT(INOUT) :: this
      INTEGER,             INTENT(IN)    :: n_maxhub !how many hubbard1 procedures can there be
      INTEGER,             INTENT(IN)    :: n_maxargs !how many additional args can there be

      this%iter = 0
      this%l_runthisiter = .false.
      this%beta = 100.0   !Default temperature 1/eV

      ALLOCATE (this%init_occ(n_maxhub))
      ALLOCATE (this%l_soc_given(n_maxhub))
      ALLOCATE (this%l_ccf_given(n_maxhub))
      this%l_soc_given = .false.
      this%l_ccf_given = .false.

      ALLOCATE (this%n_addArgs(n_maxhub))
      ALLOCATE (this%arg_keys(n_maxhub,n_maxargs))
      ALLOCATE (this%arg_vals(n_maxhub,n_maxargs))
      this%n_addArgs = 0

      ALLOCATE (this%n_exc_given(n_maxhub))
      ALLOCATE (this%exc(n_maxhub,lmaxU_const-1))
      ALLOCATE (this%exc_l(n_maxhub,lmaxU_const-1))
      ALLOCATE (this%mag_mom(n_maxhub,lmaxU_const-1))
      ALLOCATE (this%init_mom(n_maxhub,lmaxU_const-1))
      this%n_exc_given = 0.0
      this%exc = 0.0
      this%exc_l = -1
      this%mag_mom = 0.0
      this%init_mom = 0.0

      ALLOCATE (this%xi(n_maxhub))
      ALLOCATE (this%ccf(n_maxhub))
      ALLOCATE (this%ccfmat(n_maxhub,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
      this%xi = 0.0
      this%ccf = 0.0
      this%ccfmat = 0.0

   END SUBROUTINE

END MODULE m_types_hub1ham