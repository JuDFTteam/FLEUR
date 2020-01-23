MODULE m_types_hybdat
  IMPLICIT NONE


   TYPE t_hybdat
      LOGICAL                ::  l_subvxc = .false.
      LOGICAL                ::  l_calhf = .false.
      LOGICAL                ::  l_addhf = .false.
      INTEGER                :: lmaxcd, maxindxc
      INTEGER                :: maxfac
      REAL, ALLOCATABLE      :: gridf(:,:)
      INTEGER, ALLOCATABLE   :: nindxc(:,:)
      INTEGER, ALLOCATABLE   :: lmaxc(:)
      REAL, ALLOCATABLE      :: core1(:,:,:,:), core2(:,:,:,:)
      REAL, ALLOCATABLE      :: eig_c(:,:,:)
      INTEGER, ALLOCATABLE   :: kveclo_eig(:,:)
      REAL, ALLOCATABLE      :: sfac(:), fac(:)
      REAL, ALLOCATABLE      :: gauntarr(:,:,:,:,:,:)
      REAL, ALLOCATABLE      :: bas1(:,:,:,:), bas2(:,:,:,:)
      REAL, ALLOCATABLE      :: bas1_MT(:,:,:), drbas1_MT(:,:,:)
      REAL, ALLOCATABLE      :: prodm(:,:,:,:)
      REAL, ALLOCATABLE      :: div_vv(:, :, :)
      INTEGER, ALLOCATABLE   :: pntgptd(:)
      INTEGER, ALLOCATABLE   :: pntgpt(:,:,:,:)
      INTEGER, ALLOCATABLE   :: nindxp1(:,:)
      INTEGER, ALLOCATABLE   ::  ne_eig(:)
      INTEGER, ALLOCATABLE   ::  nbands(:)
      INTEGER, ALLOCATABLE   ::  nobd(:, :)
      INTEGER                ::  maxlmindx = -1
      COMPLEX, ALLOCATABLE   :: stepfunc(:,:,:)
      INTEGER                ::  nbasp = -1
      INTEGER                ::  maxbasm1 = -1
      INTEGER, ALLOCATABLE   ::  nbasm(:)
   contains
      procedure :: set_stepfunction => set_stepfunction
      procedure :: free => free_hybdat
      procedure :: allocate => allocate_hybdat
   END TYPE t_hybdat

contains
   subroutine allocate_hybdat(hybdat, atoms, num_radfun_per_l)
      use m_types_atoms
      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_atoms), intent(in)      :: atoms
      integer, intent(in)            :: num_radfun_per_l(:,:)

      allocate(hybdat%lmaxc(atoms%ntype), source=0)
      allocate(hybdat%bas1(atoms%jmtd, maxval(num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)
      allocate(hybdat%bas2(atoms%jmtd, maxval(num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)
      allocate(hybdat%bas1_MT(maxval(num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)
      allocate(hybdat%drbas1_MT(maxval(num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)
   end subroutine allocate_hybdat

   subroutine free_hybdat(hybdat)
      implicit none
      class(t_hybdat), intent(inout) :: hybdat

      if(allocated(hybdat%lmaxc)) deallocate(hybdat%lmaxc)
      if(allocated(hybdat%bas1)) deallocate(hybdat%bas1)
      if(allocated(hybdat%bas2)) deallocate(hybdat%bas2)
      if(allocated(hybdat%bas1_MT)) deallocate(hybdat%bas1_MT)
      if(allocated(hybdat%drbas1_MT)) deallocate(hybdat%drbas1_MT)
      if(allocated(hybdat%nindxc)) deallocate(hybdat%nindxc)
      if(allocated(hybdat%core1)) deallocate(hybdat%core1)
      if(allocated(hybdat%core2)) deallocate(hybdat%core2)
      if(allocated(hybdat%eig_c)) deallocate(hybdat%eig_c)
      if(allocated(hybdat%fac)) deallocate(hybdat%fac)
      if(allocated(hybdat%gauntarr)) deallocate(hybdat%gauntarr)
   end subroutine free_hybdat

   subroutine set_stepfunction(hybdat, cell, atoms, g, svol)
      use m_types_cell
      use m_types_atoms
      use m_judft
      implicit none
      class(t_hybdat),INTENT(INOUT) :: hybdat
      type(t_cell),  INTENT(in)    :: cell
      type(t_atoms), INTENT(in)    :: atoms
      integer,       INTENT(in)    :: g(3)
      real,          INTENT(in)    :: svol
      integer :: i, j, k, ok

      if (.not. allocated(hybdat%stepfunc)) then
         call timestart("setup stepfunction")
         ALLOCATE (hybdat%stepfunc(-g(1):g(1), -g(2):g(2), -g(3):g(3)), stat=ok)
         IF (ok /= 0) then
            call juDFT_error('wavefproducts_inv5: error allocation stepfunc')
         endif

         DO i = -g(1), g(1)
            DO j = -g(2), g(2)
               DO k = -g(3), g(3)
                  hybdat%stepfunc(i,j,k) = stepfunction(cell, atoms, [i, j, k])/svol
               END DO
            END DO
         END DO
         call timestop("setup stepfunction")
      endif

   end subroutine set_stepfunction

   !private subroutine
   FUNCTION stepfunction(cell, atoms, g)
      USE m_types_cell
      USE m_types_atoms
      USE m_constants
      IMPLICIT NONE

      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_atoms), INTENT(IN)   :: atoms

      INTEGER, INTENT(IN) :: g(3)
      COMPLEX             :: stepfunction  !Is real in inversion case
      REAL                :: gnorm, gnorm3, r, fgr
      INTEGER             :: itype, ieq, icent

      gnorm = gptnorm(g, cell%bmat)
      gnorm3 = gnorm**3
      IF (abs(gnorm) < 1e-12) THEN
         stepfunction = 1
         DO itype = 1, atoms%ntype
            stepfunction = stepfunction - atoms%neq(itype)*atoms%volmts(itype)/cell%omtil
         END DO
      ELSE
         stepfunction = 0
         icent = 0
         DO itype = 1, atoms%ntype
            r = gnorm*atoms%rmt(itype)
            fgr = fpi_const*(sin(r) - r*cos(r))/gnorm3/cell%omtil
            DO ieq = 1, atoms%neq(itype)
               icent = icent + 1
               stepfunction = stepfunction - fgr*exp(-cmplx(0., tpi_const*dot_product(atoms%taual(:,icent), g)))
            ENDDO
         ENDDO
      ENDIF

   END FUNCTION stepfunction

   PURE FUNCTION gptnorm(gpt, bmat)
      IMPLICIT NONE
      REAL                :: gptnorm
      INTEGER, INTENT(IN)  :: gpt(3)
      REAL, INTENT(IN)     :: bmat(3, 3)

      gptnorm = norm2(matmul(gpt(:), bmat(:,:)))

   END FUNCTION gptnorm

 END MODULE m_types_hybdat
