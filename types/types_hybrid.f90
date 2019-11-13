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
      INTEGER                ::  maxbasm1
      INTEGER                ::  max_indx_p_1
      INTEGER                ::  maxlmindx
      INTEGER, ALLOCATABLE   ::  num_radfun_per_l(:,:)
      INTEGER, ALLOCATABLE   ::  select1(:,:)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:,:)
      INTEGER, ALLOCATABLE   ::  tvec(:,:,:)
      INTEGER, ALLOCATABLE   ::  nbasm(:)
      !REAL, ALLOCATABLE      ::  radbasfn_mt(:,:,:,:)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:,:,:,:)
      INTEGER, ALLOCATABLE   ::  ne_eig(:)
      INTEGER, ALLOCATABLE   ::  nbands(:)
      INTEGER, ALLOCATABLE   ::  nobd(:,:)
      REAL, ALLOCATABLE      ::  div_vv(:,:,:)
   CONTAINS
      procedure :: set_num_radfun_per_l => set_num_radfun_per_l_hybrid
   END TYPE t_hybrid

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
      INTEGER, ALLOCATABLE   :: pntgptd(:)
      INTEGER, ALLOCATABLE   :: pntgpt(:,:,:,:)
      INTEGER, ALLOCATABLE   :: nindxp1(:,:)
      COMPLEX, ALLOCATABLE   :: stepfunc(:,:,:)
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
