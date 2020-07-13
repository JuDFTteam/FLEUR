MODULE m_types_hybdat
   use m_types_usdus
   use m_types_mat
   IMPLICIT NONE

   type t_coul
      REAL, ALLOCATABLE      :: mt1_r(:, :, :, :)
      REAL, ALLOCATABLE      :: mt2_r(:, :, :, :),   mt3_r(:, :, :)
      REAL, ALLOCATABLE      :: mtir_r(:, :)
      COMPLEX, ALLOCATABLE   :: mt1_c(:, :, :, :)
      COMPLEX, ALLOCATABLE   :: mt2_c(:, :, :, :),   mt3_c(:, :, :)
      COMPLEX, ALLOCATABLE   :: mtir_c(:, :)
      integer  :: bcast_req(4)
      logical  :: bcast_finished = .False.
   contains 
      procedure :: init     => t_coul_init
      procedure :: alloc    => t_coul_alloc
      procedure :: free     => t_coul_free
      procedure :: mpi_ibc  => t_coul_mpi_ibc
      procedure :: mpi_wait => t_coul_mpi_wait
      procedure :: size_MB  => t_coul_size_MB
   end type t_coul

   TYPE t_hybdat
      LOGICAL                :: l_subvxc = .false.
      LOGICAL                :: l_calhf = .false.
      LOGICAL                :: l_addhf = .false.
      LOGICAL                :: l_print_iob_splitting = .True.
      INTEGER                :: lmaxcd, maxindxc
      INTEGER                :: maxfac
      REAL, ALLOCATABLE      :: gridf(:, :)
      INTEGER, ALLOCATABLE   :: nindxc(:, :)
      INTEGER, ALLOCATABLE   :: lmaxc(:)
      REAL, ALLOCATABLE      :: core1(:, :, :, :), core2(:, :, :, :)
      REAL, ALLOCATABLE      :: eig_c(:, :, :)
      INTEGER, ALLOCATABLE   :: kveclo_eig(:, :)
      REAL, ALLOCATABLE      :: sfac(:), fac(:)
      REAL, ALLOCATABLE      :: gauntarr(:, :, :, :, :, :)
      REAL, ALLOCATABLE      :: bas1(:, :, :, :), bas2(:, :, :, :)
      REAL, ALLOCATABLE      :: bas1_MT(:, :, :), drbas1_MT(:, :, :)
      REAL, ALLOCATABLE      :: prodm(:, :, :, :)
      REAL, ALLOCATABLE      :: div_vv(:, :, :)
      INTEGER, ALLOCATABLE   :: pntgptd(:)
      INTEGER, ALLOCATABLE   :: pntgpt(:, :, :, :)
      INTEGER, ALLOCATABLE   :: nindxp1(:, :)
      INTEGER, ALLOCATABLE   :: ne_eig(:)
      INTEGER, ALLOCATABLE   :: nbands(:)
      INTEGER, ALLOCATABLE   :: nobd(:, :)
      INTEGER                :: maxlmindx = -1
      COMPLEX, ALLOCATABLE   :: stepfunc(:, :, :)
      INTEGER                :: nbasp = -1
      INTEGER                :: eig_id = -1
      INTEGER, ALLOCATABLE   :: nbasm(:)

      ! coulomb matrix stuff
      type(t_coul), allocatable :: coul(:)

      type(t_usdus)             :: usdus
      type(t_mat), allocatable  :: v_x(:,:) ! nkpt, jsp
   contains
      procedure :: set_stepfunction => set_stepfunction
      procedure :: free       => free_hybdat
      procedure :: allocate   => allocate_hybdat
      procedure :: set_states => set_states_hybdat
   END TYPE t_hybdat

contains
   subroutine set_states_hybdat(hybdat, fi, results, jsp)
      use m_judft
      use m_types_misc
      use m_types_fleurinput
      implicit none 
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in) :: fi
      type(t_results), intent(in)    :: results
      integer, intent(in)            :: jsp

      integer :: nk, i

      DO nk = 1, fi%kpts%nkpt, 1
         hybdat%ne_eig(nk) = results%neig(nk, jsp)
         hybdat%nobd(nk,jsp) = COUNT(results%w_iks(:hybdat%ne_eig(nk), nk, jsp) > 0.0)
         hybdat%nbands(nk) = fi%hybinp%bands1
      END do

      DO nk = 1, fi%kpts%nkpt
         IF (hybdat%nobd(nk,jsp) > hybdat%nbands(nk)) THEN
            hybdat%nbands(nk) = hybdat%nobd(nk,jsp)
         END IF
      ENDDO
      
      ! spread hybdat%nobd from IBZ to whole BZ
      DO nk = 1, fi%kpts%nkptf
         i = fi%kpts%bkp(nk)
         hybdat%nbands(nk) = hybdat%nbands(i)
         hybdat%nobd(nk,jsp) = hybdat%nobd(i,jsp)
      END DO
   end subroutine set_states_hybdat

   function t_coul_size_MB(coul) result(size_MB)
      implicit none 
      class(t_coul), intent(in) :: coul
      real  :: size_MB 

      size_MB = 0
      
      ! real parts
      if(allocated(coul%mt1_r))   size_MB = size_MB + 8 * 1e-6 * size(coul%mt1_r) 
      if(allocated(coul%mt2_r))   size_MB = size_MB + 8 * 1e-6 * size(coul%mt2_r) 
      if(allocated(coul%mt3_r))   size_MB = size_MB + 8 * 1e-6 * size(coul%mt3_r) 
      if(allocated(coul%mtir_r))  size_MB = size_MB + 8 * 1e-6 * size(coul%mtir_r) 

      ! complex parts
      if(allocated(coul%mt1_c))   size_MB = size_MB + 16 * 1e-6 * size(coul%mt1_c) 
      if(allocated(coul%mt2_c))   size_MB = size_MB + 16 * 1e-6 * size(coul%mt2_r) 
      if(allocated(coul%mt3_c))   size_MB = size_MB + 16 * 1e-6 * size(coul%mt3_r) 
      if(allocated(coul%mtir_c))  size_MB = size_MB + 16 * 1e-6 * size(coul%mtir_r) 
   end function

   subroutine t_coul_mpi_wait(coul)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none 
      class(t_coul)                 :: coul
#ifdef CPP_MPI
      integer :: ierr, i 
      
      if(.not. coul%bcast_finished) then      
         do i = 1,size(coul%bcast_req)
            call MPI_WAIT(coul%bcast_req(i), MPI_STATUS_IGNORE, ierr)
            if(ierr /= 0) call judft_error("Error in MPI_wait for coul%bcast_req no. " // int2str(i))
         enddo
         coul%bcast_finished = .True.
      endif 
#else
      coul%bcast_finished = .True.
#endif
   end subroutine t_coul_mpi_wait

   subroutine t_coul_mpi_ibc(coul, fi, glob_mpi, root)
      use m_types_fleurinput
      use m_types_hybmpi
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none 
      class(t_coul)                  :: coul
      type(t_fleurinput), intent(in) :: fi
      type(t_hybmpi), intent(in)     :: glob_mpi
      integer, intent(in)            :: root
#ifdef CPP_MPI
      integer :: ierr



      if (fi%sym%invs) THEN
         call MPI_IBcast(coul%mt1_r, size(coul%mt1_r), MPI_DOUBLE_PRECISION, root, glob_mpi%comm, coul%bcast_req(1), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mt1_r failed")

         call MPI_IBcast(coul%mt2_r,   size(coul%mt2_r),   MPI_DOUBLE_PRECISION, root, glob_mpi%comm, coul%bcast_req(2), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mt2_r failed")

         call MPI_IBcast(coul%mt3_r,   size(coul%mt3_r),   MPI_DOUBLE_PRECISION, root, glob_mpi%comm, coul%bcast_req(3), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mt3_r failed")

         call MPI_IBcast(coul%mtir_r,  size(coul%mtir_r),  MPI_DOUBLE_PRECISION, root, glob_mpi%comm, coul%bcast_req(4), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mtir_r failed")
      else 
         call MPI_IBcast(coul%mt1_c, size(coul%mt1_c),     MPI_DOUBLE_COMPLEX,  root, glob_mpi%comm, coul%bcast_req(1), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mt1_c failed")

         call MPI_IBcast(coul%mt2_c,   size(coul%mt2_c),   MPI_DOUBLE_COMPLEX , root, glob_mpi%comm, coul%bcast_req(2), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mt2_r failed")

         call MPI_IBcast(coul%mt3_c,   size(coul%mt3_c),   MPI_DOUBLE_COMPLEX , root, glob_mpi%comm, coul%bcast_req(3), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mt3_r failed")

         call MPI_IBcast(coul%mtir_c,  size(coul%mtir_c),  MPI_DOUBLE_COMPLEX , root, glob_mpi%comm, coul%bcast_req(4), ierr)
         if(ierr /= 0) call judft_error("MPI_IBcast of coul%mtir_r failed")
      endif
#endif
      coul%bcast_finished = .False.
   end subroutine t_coul_mpi_ibc

   subroutine t_coul_free(coul)
      implicit none 
      class(t_coul), intent(inout) :: coul 

      if(allocated(coul%mt1_r)) deallocate(coul%mt1_r)
      if(allocated(coul%mt1_c)) deallocate(coul%mt1_c)
      if(allocated(coul%mt2_r)) deallocate(coul%mt2_r)
      if(allocated(coul%mt3_r)) deallocate(coul%mt3_r)
      if(allocated(coul%mtir_r)) deallocate(coul%mtir_r)
      if(allocated(coul%mt2_c)) deallocate(coul%mt2_c)
      if(allocated(coul%mt3_c)) deallocate(coul%mt3_c)
      if(allocated(coul%mtir_c)) deallocate(coul%mtir_c)
   end subroutine t_coul_free

   subroutine t_coul_alloc(coul, fi, num_radbasfn, n_g, ikpt)
      use m_types_fleurinput
      use m_judft
      implicit NONE 
      class(t_coul), intent(inout) :: coul
      type(t_fleurinput), intent(in)    :: fi
      integer, intent(in) :: num_radbasfn(:, :), n_g(:), ikpt
      integer :: info, isize, l, itype

      isize = sum([(((2*l + 1)*fi%atoms%neq(itype), l=0, fi%hybinp%lcutm1(itype)),&
                                                    itype=1, fi%atoms%ntype)]) &
                  + n_g(ikpt)



      if (fi%sym%invs) THEN      
         if(.not. allocated(coul%mt1_r)) then 
            allocate(coul%mt1_r(maxval(num_radbasfn) - 1,&
                              maxval(num_radbasfn) - 1,&
                              0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt1_r")
         endif
         if(.not. allocated(coul%mt2_r)) then
            allocate(coul%mt2_r(maxval(num_radbasfn) - 1,&
                                 -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1),&
                                 0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt2_r")
         endif

         if(.not. allocated(coul%mt3_r)) then
             allocate(coul%mt3_r(maxval(num_radbasfn) - 1,fi%atoms%nat, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt3_r")
         endif 

         if(.not. allocated(coul%mtir_r)) then
            allocate(coul%mtir_r(isize,isize), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mtir_r")
         endif 
      else      
         if(.not. allocated(coul%mt1_c)) then 
            allocate(coul%mt1_c(maxval(num_radbasfn) - 1,&
                              maxval(num_radbasfn) - 1,&
                              0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt1_c")
         endif
         if(.not. allocated(coul%mt2_c)) then
            allocate(coul%mt2_c(maxval(num_radbasfn) - 1,&
                                 -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), &
                                 0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt2_c")
         endif 

         if(.not. allocated(coul%mt3_c)) then
            allocate(coul%mt3_c(maxval(num_radbasfn) - 1, fi%atoms%nat, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt3_c")
         endif 

         if(.not. allocated(coul%mtir_c)) then
            allocate(coul%mtir_c(isize, isize), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mtir_c")
         endif 
      endif
   end subroutine t_coul_alloc

   subroutine t_coul_init(coul)
      implicit none 
      class(t_coul), intent(inout) :: coul
      
      if(allocated(coul%mt1_r)) coul%mt1_r = 0
      if(allocated(coul%mt1_c)) coul%mt1_c = 0

      if(allocated(coul%mt2_r))   coul%mt2_r = 0
      if(allocated(coul%mt3_r))   coul%mt3_r = 0
      
      if(allocated(coul%mt2_c))   coul%mt2_c = 0
      if(allocated(coul%mt3_c))   coul%mt3_c = 0
   end subroutine t_coul_init

   subroutine allocate_hybdat(hybdat, fi, num_radfun_per_l)
      use m_types_fleurinput
      use m_judft
      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: num_radfun_per_l(:, :)
      integer                        :: ok(12)

      ok = -1
      allocate(hybdat%lmaxc(fi%atoms%ntype), &
               stat=ok(1), source=0)
      allocate(hybdat%bas1(fi%atoms%jmtd, maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(2), source=0.0)
      allocate(hybdat%bas2(fi%atoms%jmtd, maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(3), source=0.0)
      allocate(hybdat%bas1_MT(maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(4), source=0.0)
      allocate(hybdat%drbas1_MT(maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(5), source=0.0)

      ! core allocs
      allocate(hybdat%nindxc(0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(6), source=0)
      allocate(hybdat%core1(fi%atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(7), source=0.0)
      allocate(hybdat%core2(fi%atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(8), source=0.0)
      allocate(hybdat%eig_c(hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(9), source=0.0)

      allocate(hybdat%fac(0:hybdat%maxfac), stat=ok(10), source=0.0)
      allocate(hybdat%sfac(0:hybdat%maxfac), stat=ok(11), source=0.0)

      ALLOCATE(hybdat%gauntarr(2, 0:fi%atoms%lmaxd, 0:fi%atoms%lmaxd, 0:maxval(fi%hybinp%lcutm1), &
                               -fi%atoms%lmaxd:fi%atoms%lmaxd, -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1)), &
               stat=ok(12), source=0.0)

      if(any(ok /= 0)) then
         write(*, *) "allocation of hybdat failed. Error in array no.:"
         write(*, *) maxloc(abs(ok))
         call juDFT_error("allocation of hybdat failed. Error in array no is the outfile")
      endif
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
      if(allocated(hybdat%sfac)) deallocate(hybdat%sfac)
      if(allocated(hybdat%gauntarr)) deallocate(hybdat%gauntarr)
   end subroutine free_hybdat

   subroutine set_stepfunction(hybdat, cell, atoms, g, svol)
      use m_types_cell
      use m_types_atoms
      use m_judft
      implicit none
      class(t_hybdat), INTENT(INOUT) :: hybdat
      type(t_cell), INTENT(in)    :: cell
      type(t_atoms), INTENT(in)    :: atoms
      integer, INTENT(in)    :: g(3)
      real, INTENT(in)    :: svol
      integer :: i, j, k, ok

      if(.not. allocated(hybdat%stepfunc)) then
         call timestart("setup stepfunction")
         ALLOCATE(hybdat%stepfunc(-g(1):g(1), -g(2):g(2), -g(3):g(3)), stat=ok)
         IF(ok /= 0) then
            call juDFT_error('wavefproducts_inv5: error allocation stepfunc')
         endif

         DO i = -g(1), g(1)
            DO j = -g(2), g(2)
               DO k = -g(3), g(3)
                  hybdat%stepfunc(i, j, k) = stepfunction(cell, atoms,[i, j, k])/svol
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
      IF(abs(gnorm) < 1e-12) THEN
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
               stepfunction = stepfunction - fgr*exp(-cmplx(0., tpi_const*dot_product(atoms%taual(:, icent), g)))
            ENDDO
         ENDDO
      ENDIF

   END FUNCTION stepfunction

   PURE FUNCTION gptnorm(gpt, bmat)
      IMPLICIT NONE
      REAL                :: gptnorm
      INTEGER, INTENT(IN)  :: gpt(3)
      REAL, INTENT(IN)     :: bmat(3, 3)

      gptnorm = norm2(matmul(gpt(:), bmat(:, :)))

   END FUNCTION gptnorm

END MODULE m_types_hybdat
