MODULE m_types_hybrid
   IMPLICIT NONE

   TYPE t_hybrid
      LOGICAL                ::  l_hybrid = .false.
      LOGICAL                ::  l_subvxc = .false.
      LOGICAL                ::  l_calhf = .false.
      LOGICAL                ::  l_addhf = .false.
      INTEGER                ::  ewaldlambda
      INTEGER                ::  lexp = 0
      INTEGER                ::  bands1 !Only read in
      INTEGER                ::  nbasp
      INTEGER                ::  maxlcutm1
      INTEGER                ::  maxbasm1
      INTEGER                ::  max_indx_p_1
      INTEGER                ::  maxgptm1
      INTEGER                ::  maxlmindx
      INTEGER                ::  gptmd
      INTEGER, ALLOCATABLE   ::  num_radfun_per_l(:,:)
      INTEGER, ALLOCATABLE   ::  select1(:,:)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  nindxm1(:,:)
      INTEGER, ALLOCATABLE   ::  gptm(:,:)
      INTEGER, ALLOCATABLE   ::  ngptm1(:)
      INTEGER, ALLOCATABLE   ::  pgptm1(:,:)
      INTEGER, ALLOCATABLE   ::  ngptm(:)
      INTEGER, ALLOCATABLE   ::  pgptm(:,:)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:,:)
      INTEGER, ALLOCATABLE   ::  tvec(:,:,:)
      INTEGER, ALLOCATABLE   ::  nbasm(:)
      REAL                   ::  gcutm
      REAL                   ::  tolerance1  !only read in
      REAL, ALLOCATABLE      ::  basm1(:,:,:,:)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:,:,:,:)
      INTEGER, ALLOCATABLE   ::  ne_eig(:), nbands(:), nobd(:)
      REAL, ALLOCATABLE      ::  div_vv(:,:,:)
   CONTAINS
      procedure :: set_num_radfun_per_l => set_num_radfun_per_l_hybrid
   END TYPE t_hybrid

   TYPE t_prodtype
      INTEGER, ALLOCATABLE :: l1(:,:,:)
      INTEGER, ALLOCATABLE :: l2(:,:,:)
      INTEGER, ALLOCATABLE :: n1(:,:,:)
      INTEGER, ALLOCATABLE :: n2(:,:,:)
   contains
      procedure  :: init   => init_prodtype
      procedure  :: free   => free_prodtype
      procedure  :: set_nl => set_nl_prodtype
   END TYPE t_prodtype

   TYPE t_hybdat
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
      TYPE(t_PRODTYPE)       :: prod
      INTEGER, ALLOCATABLE   :: pntgptd(:)
      INTEGER, ALLOCATABLE   :: pntgpt(:,:,:,:)
      INTEGER, ALLOCATABLE   :: nindxp1(:,:)
      COMPLEX, ALLOCATABLE   ::  stepfunc(:,:,:)
   contains
      procedure  :: set_stepfunction => set_stepfunction
   END TYPE t_hybdat

contains
   subroutine set_stepfunction(hybdat, cell, atoms, g, svol)
      use m_types_setup
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
      USE m_types_setup
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

   subroutine init_prodtype(prod, hybrid, atoms)
      use m_types_setup
      use m_judft
      implicit none
      class(t_prodtype)          :: prod
      type(t_hybrid), intent(in) :: hybrid
      type(t_atoms),  intent(in) :: atoms
      integer                    :: ok

      ALLOCATE (prod%l1(hybrid%max_indx_p_1, 0:hybrid%maxlcutm1, atoms%ntype), stat=ok)
      IF (ok /= 0) call judft_error('init_prodtype: failure allocation prod%l1')

      ALLOCATE (prod%l2, mold=prod%l1, stat=ok)
      IF (ok /= 0) call judft_error('init_prodtype: failure allocation prod%l2')

      ALLOCATE (prod%n1, mold=prod%l1, stat=ok)
      IF (ok /= 0) call judft_error('init_prodtype: failure allocation prod%n1')

      ALLOCATE (prod%n2, mold=prod%l1, stat=ok)
      IF (ok /= 0) call judft_error('init_prodtype: failure allocation prod%n2')
   end subroutine init_prodtype

   subroutine free_prodtype(prod)
      use m_types_setup
      implicit NONE
      class(t_prodtype)          :: prod

      IF(ALLOCATED(prod%l1)) DEALLOCATE (prod%l1)
      IF(ALLOCATED(prod%l2)) DEALLOCATE (prod%l2)
      IF(ALLOCATED(prod%n1)) DEALLOCATE (prod%n1)
      IF(ALLOCATED(prod%n2)) DEALLOCATE (prod%n2)
   end subroutine free_prodtype

   subroutine set_nl_prodtype(prod,n,l,itype,n1,l1,n2,l2)
      use m_types_setup
      implicit NONE
      class(t_prodtype)    :: prod
      integer, intent(in)  :: n, l, itype
      integer, intent(out) :: n1, l1, n2, l2

      l1 = prod%l1(n,l,itype) !
      l2 = prod%l2(n,l,itype) ! current basis-function product
      n1 = prod%n1(n,l,itype) ! = bas(:,n1,l1,itype)*bas(:,n2,l2,itype) = b1*b2
      n2 = prod%n2(n,l,itype) !

   end subroutine set_nl_prodtype

   subroutine set_num_radfun_per_l_hybrid(hybrid, atoms)
      use m_types_setup
      implicit NONE
      class(t_hybrid) :: hybrid
      type(t_atoms)   :: atoms
      integer :: itype, ilo

      ! there is always at least two: u and u_dot
      hybrid%num_radfun_per_l = 2
      DO itype = 1, atoms%ntype
         DO ilo = 1, atoms%nlo(itype)
            hybrid%num_radfun_per_l(atoms%llo(ilo, itype), itype) &
              = hybrid%num_radfun_per_l(atoms%llo(ilo, itype), itype) + 1
         END DO
      END DO
   end subroutine set_num_radfun_per_l_hybrid
END MODULE m_types_hybrid
