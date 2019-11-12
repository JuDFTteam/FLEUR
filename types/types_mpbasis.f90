module m_types_mpbasis
   implicit none

   type t_mpbasis
      INTEGER, ALLOCATABLE   :: gptm(:,:) ! (3, num_gpts)
      INTEGER, ALLOCATABLE   :: ngptm(:)
      INTEGER, ALLOCATABLE   :: gptm_ptr(:,:)
      REAL                   :: g_cutoff
      INTEGER, ALLOCATABLE   :: num_radbasfn(:,:)
      REAL, ALLOCATABLE      :: radbasfn_mt(:,:,:,:)
      REAL                   :: linear_dep_tol  !only read in
   CONTAINS
      procedure :: num_gpts               => mpbasis_num_gpts
      procedure :: gen_gvec               => mpbasis_gen_gvec
      procedure :: check_orthonormality   => mpbasis_check_orthonormality
      procedure :: check_radbasfn         => mpbasis_check_radbasfn
      procedure :: calc_olap_radbasfn     => mpbasis_calc_olap_radbasfn
      procedure :: filter_radbasfn        => mpbasis_filter_radbasfn
      procedure :: trafo_to_orthonorm_bas => mpbasis_trafo_to_orthonorm_bas
      procedure :: add_l0_fun             => mpbasis_add_l0_fun
      procedure :: reduce_linear_dep      => mpbasis_reduce_linear_dep
   end type t_mpbasis
contains
   function mpbasis_num_gpts(mpbasis)
      implicit NONE
      class(t_mpbasis), intent(in) :: mpbasis
      integer    :: mpbasis_num_gpts

      mpbasis_num_gpts = size(mpbasis%gptm,dim=2)
   end function mpbasis_num_gpts

   subroutine mpbasis_gen_gvec(mpbasis, cell, kpts, mpi)
      use m_types_setup
      use m_types_kpts
      use m_types_mpi
      USE m_intgrf, ONLY: intgrf_init, intgrf
      use m_rorder, only: rorderpf
      implicit NONE
      class(t_mpbasis), intent(inout) :: mpbasis
      type(t_cell), intent(in)       :: cell
      type(t_kpts), intent(in)       :: kpts
      type(t_mpi), intent(in)        :: mpi


      integer :: i, n, n1, n2, divconq
      integer :: x, y, z, ikpt, igpt
      integer :: g(3)
      real    :: longest_k, kvec(3)
      logical :: l_found_new_gpt, l_found_kg_in_sphere

      INTEGER, ALLOCATABLE            ::  unsrt_pgptm(:,:) ! unsorted pointers to g vectors
      REAL, ALLOCATABLE               ::  length_kg(:,:) ! length of the vectors k + G
      INTEGER, ALLOCATABLE            ::  ptr(:)

      allocate(mpbasis%ngptm(kpts%nkptf))

      mpbasis%ngptm = 0
      i = 0
      n = -1

      longest_k = MAXVAL([(norm2(MATMUL(kpts%bkf(:,ikpt), cell%bmat)), ikpt=1, kpts%nkptf)])

      ! a first run for the determination of the dimensions of the fields gptm,pgptm

      do
         n = n + 1
         l_found_new_gpt = .FALSE.
         do x = -n, n
            n1 = n - ABS(x)
            do y = -n1, n1
               n2 = n1 - ABS(y)
               do z = -n2, n2, MAX(2*n2, 1)
                  g = [x, y, z]
                  IF ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpbasis%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  do ikpt = 1, kpts%nkptf
                     IF (norm2(MATMUL(kpts%bkf(:,ikpt) + g, cell%bmat)) <= mpbasis%g_cutoff) THEN
                        IF (.NOT. l_found_kg_in_sphere) THEN
                           i = i + 1
                           l_found_kg_in_sphere = .TRUE.
                        END IF

                        mpbasis%ngptm(ikpt) = mpbasis%ngptm(ikpt) + 1
                        l_found_new_gpt = .TRUE.
                     END IF
                  enddo ! k-loop
               enddo
            enddo
         enddo
         IF (.NOT. l_found_new_gpt) EXIT
      enddo

      allocate(mpbasis%gptm(3,i)) ! i = gptmd
      allocate(mpbasis%gptm_ptr(maxval(mpbasis%ngptm), kpts%nkptf))

      ! Allocate and initialize arrays needed for G vector ordering
      allocate(unsrt_pgptm(maxval(mpbasis%ngptm), kpts%nkptf))
      allocate(length_kG(maxval(mpbasis%ngptm), kpts%nkptf))

      mpbasis%gptm = 0
      mpbasis%gptm_ptr = 0
      mpbasis%ngptm = 0

      i = 0
      n = -1

      length_kG = 0
      unsrt_pgptm = 0

      do
         n = n + 1
         l_found_new_gpt = .FALSE.
         do x = -n, n
            n1 = n - ABS(x)
            do y = -n1, n1
               n2 = n1 - ABS(y)
               do z = -n2, n2, MAX(2*n2, 1)
                  g = [x, y, z]
                  IF ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpbasis%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  do ikpt = 1, kpts%nkptf
                     kvec = kpts%bkf(:,ikpt)

                     IF (norm2(MATMUL(kvec + g, cell%bmat)) <= mpbasis%g_cutoff) THEN
                        IF (.NOT. l_found_kg_in_sphere) THEN
                           i = i + 1
                           mpbasis%gptm(:,i) = g
                           l_found_kg_in_sphere = .TRUE.
                        END IF

                        mpbasis%ngptm(ikpt) = mpbasis%ngptm(ikpt) + 1
                        l_found_new_gpt = .TRUE.

                        ! Position of the vector is saved as pointer
                        unsrt_pgptm(mpbasis%ngptm(ikpt), ikpt) = i
                        ! Save length of vector k + G for array sorting
                        length_kG(mpbasis%ngptm(ikpt), ikpt) = norm2(MATMUL(kvec + g, cell%bmat))
                     END IF
                  enddo
               enddo
            enddo
         enddo
         IF (.NOT. l_found_new_gpt) EXIT
      enddo

      ! Sort pointers in array, so that shortest |k+G| comes first
      do ikpt = 1, kpts%nkptf
         allocate(ptr(mpbasis%ngptm(ikpt)))
         ! Divide and conquer algorithm for arrays > 1000 entries
         divconq = MAX(0, INT(1.443*LOG(0.001*mpbasis%ngptm(ikpt))))
         ! create pointers which correspond to a sorted array
         CALL rorderpf(ptr, length_kG(1:mpbasis%ngptm(ikpt), ikpt), mpbasis%ngptm(ikpt), divconq)
         ! rearrange old pointers
         do igpt = 1, mpbasis%ngptm(ikpt)
            mpbasis%gptm_ptr(igpt, ikpt) = unsrt_pgptm(ptr(igpt), ikpt)
         enddo
         deallocate(ptr)
      enddo

      IF (mpi%irank == 0) THEN
         WRITE (6, '(/A)') 'Mixed basis'
         WRITE (6, '(A,I5)') 'Number of unique G-vectors: ', mpbasis%num_gpts()
         WRITE (6, *)
         WRITE (6, '(3x,A)') 'IR Plane-wave basis with cutoff of gcutm (mpbasis%g_cutoff/2*input%rkmax):'
         WRITE (6, '(5x,A,I5)') 'Maximal number of G-vectors:', maxval(mpbasis%ngptm)
      END IF
   end subroutine mpbasis_gen_gvec

   subroutine mpbasis_check_orthonormality(mpbasis, atoms, mpi, l, itype, gridf)
      USE m_intgrf, ONLY: intgrf
      use m_types_setup
      use m_types_mpi
      use m_judft
      implicit none
      class(t_mpbasis)          :: mpbasis
      type(t_atoms), intent(in) :: atoms
      type(t_mpi), intent(in)   :: mpi
      integer, intent(in)       :: itype, l
      real, intent(in)          :: gridf(:,:)

      real              :: overlap, error, cum_err_sq
      integer           :: i, j, n_radbasfn

      CHARACTER, PARAMETER :: lchar(0:38) = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', &
                                             'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', &
                                             'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x']


      n_radbasfn = mpbasis%num_radbasfn(l, itype)

      cum_err_sq = 0
      do i = 1, n_radbasfn
         do j = 1, i
            overlap = intgrf(mpbasis%radbasfn_mt(:,i, l, itype)*mpbasis%radbasfn_mt(:,j, l, itype), &
                           atoms, itype, gridf)

            IF (i == j) THEN
               error = ABS(1 - overlap)
               cum_err_sq = cum_err_sq + overlap**2
            ELSE
               error = ABS(overlap)
               cum_err_sq = cum_err_sq + 2*error**2
            END IF

            IF (error > 1e-6) THEN
               IF (mpi%irank == 0) THEN
                  WRITE (6, '(A)') 'mixedbasis: Bad orthonormality of ' &
                     //lchar(l)//'-product basis. Increase tolerance.'
                  WRITE (6, '(12X,A,F9.6,A,2(I3.3,A))') 'Deviation of', &
                     error, ' occurred for (', i, ',', j, ')-overlap.'
               END IF
               CALL judft_error("Bad orthonormality of product basis", &
                                 hint='Increase tolerance', &
                                 calledby='mixedbasis%check_orthonormality')
            END IF

         enddo
      enddo

      IF (mpi%irank == 0) THEN
         WRITE (6, '(6X,A,I4,''   ('',ES8.1,'' )'')') &
                    lchar(l) // ':', n_radbasfn, SQRT(cum_err_sq)/n_radbasfn
      END IF
   end subroutine mpbasis_check_orthonormality

   subroutine mpbasis_check_radbasfn(mpbasis, atoms, hybrid)
      use m_judft
      use m_types_hybrid
      use m_types_setup
      implicit none
      class(t_mpbasis), intent(in) :: mpbasis
      type(t_atoms), intent(in)    :: atoms
      type(t_hybrid), intent(in)   :: hybrid

      integer :: itype

      do itype = 1, atoms%ntype
         IF (ANY(mpbasis%num_radbasfn(0:hybrid%lcutm1(itype), itype) == 0)) THEN
            call judft_error('any mpbasis%num_radbasfn eq 0', calledby='mixedbasis')
         endif
      enddo
   end subroutine mpbasis_check_radbasfn

   subroutine mpbasis_calc_olap_radbasfn(mpbasis, atoms, l, itype, gridf, olap)
      USE m_intgrf, ONLY: intgrf
      use m_types_setup
      implicit NONE
      class(t_mpbasis), intent(in)       :: mpbasis
      type(t_atoms), intent(in)          :: atoms
      integer, intent(in)                :: l, itype
      real, intent(in)                   :: gridf(:,:)
      real, intent(inout), allocatable   :: olap(:,:)

      integer  :: n1, n2, n_radbasfn

      n_radbasfn = mpbasis%num_radbasfn(l, itype)
      if(allocated(olap)) then
         if(any(shape(olap) /= n_radbasfn)) then
            deallocate(olap)
         endif
      endif
      if(.not. allocated(olap)) allocate(olap(n_radbasfn, n_radbasfn), source=0.0)

      DO n2 = 1, n_radbasfn
         DO n1 = 1, n2
            olap(n1, n2) = intgrf(mpbasis%radbasfn_mt(:,n1, l, itype)*mpbasis%radbasfn_mt(:,n2, l, itype), &
                                  atoms, itype, gridf)
            olap(n2, n1) = olap(n1, n2)
         END DO
      END DO
   end subroutine mpbasis_calc_olap_radbasfn

   subroutine mpbasis_filter_radbasfn(mpbasis, l, itype, n_radbasfn, eig, eigv)
      ! Get rid of linear dependencies (eigenvalue <= mpbasis%linear_dep_tol)
      implicit none
      class(t_mpbasis), intent(inout)       :: mpbasis
      integer, intent(in)                   :: l, itype, n_radbasfn
      real, intent(inout)                   :: eig(:), eigv(:,:)

      integer              :: num_radbasfn, i_bas
      integer, allocatable :: remaining_basfn(:)

      allocate(remaining_basfn(n_radbasfn), source=1)
      num_radbasfn = 0

      DO i_bas = 1, mpbasis%num_radbasfn(l, itype)
         IF (eig(i_bas) > mpbasis%linear_dep_tol) THEN
            num_radbasfn = num_radbasfn + 1
            remaining_basfn(num_radbasfn) = i_bas
         END IF
      END DO

      mpbasis%num_radbasfn(l, itype) = num_radbasfn
      eig = eig(remaining_basfn)
      eigv(:,:) = eigv(:,remaining_basfn)
   end subroutine mpbasis_filter_radbasfn

   subroutine mpbasis_diagonialize_olap(olap, eig_val, eig_vec)
      use m_judft
      implicit NONE
      real, intent(in)  :: olap(:,:)
      real, allocatable :: eig_val(:), eig_vec(:,:)

      integer              :: n, size_iwork, info
      real                 :: size_work
      integer, allocatable :: iwork(:)
      real, allocatable    :: work(:)

      if(size(olap, dim=1) /= size(olap, dim=2)) then
         call juDFT_error("only square matrices can be diagonalized")
      endif

      n = size(olap, dim=1)

      if(size(eig_val) /= n) deallocate(eig_val)
      if(.not. allocated(eig_val)) allocate(eig_val(n))

      eig_vec = olap
      ! get sizes of work arrays
      call dsyevd('V', 'U', n, eig_vec, n, eig_val,&
                  size_work, -1, size_iwork, -1, info)
      if(info /= 0) call juDFT_error("diagonalization for size failed")

      allocate(work(int(size_work)))
      allocate(iwork(size_iwork))

      call dsyevd('V', 'U', n, eig_vec, n, eig_val,&
                  work, int(size_work), iwork, size_iwork, info)
      if(info /= 0) call juDFT_error("diagonalization failed")

   end subroutine mpbasis_diagonialize_olap

   subroutine mpbasis_trafo_to_orthonorm_bas(mpbasis, full_n_radbasfn, n_grid_pt, l, itype, eig, eigv)
      implicit NONE
      class(t_mpbasis), intent(inout)  :: mpbasis
      integer, intent(in)              :: full_n_radbasfn, n_grid_pt, l, itype
      real, intent(in)                 :: eig(:), eigv(:,:)

      integer :: nn, i

      ! reduced number of basis functions
      nn = mpbasis%num_radbasfn(l, itype)

      DO i = 1, n_grid_pt
         mpbasis%radbasfn_mt(i, 1:nn, l, itype) &
            = MATMUL(mpbasis%radbasfn_mt(i, 1:full_n_radbasfn, l, itype), eigv(:,1:nn))/SQRT(eig(:nn))
      END DO
   end subroutine mpbasis_trafo_to_orthonorm_bas

   subroutine mpbasis_add_l0_fun(mpbasis, atoms, hybrid, n_grid_pt, l, itype, gridf)
      use m_types_setup
      use m_types_hybrid
      USE m_intgrf, ONLY: intgrf
      implicit none
      class(t_mpbasis), intent(inout) :: mpbasis
      type(t_atoms), intent(in)       :: atoms
      type(t_hybrid), intent(in)      :: hybrid
      integer, intent(in)             :: n_grid_pt, l, itype
      real, intent(in)                :: gridf(:,:)

      integer                         :: i, j, nn
      REAL, ALLOCATABLE               :: basmhlp(:,:,:,:)
      real                            :: norm

      nn = mpbasis%num_radbasfn(l, itype)
      IF (l == 0) THEN

         ! Check if radbasfn_mt must be reallocated
         IF (nn + 1 > SIZE(mpbasis%radbasfn_mt, 2)) THEN
            allocate(basmhlp(atoms%jmtd, nn + 1, 0:maxval(hybrid%lcutm1), atoms%ntype))
            basmhlp(:,1:nn, :,:) = mpbasis%radbasfn_mt
            deallocate(mpbasis%radbasfn_mt)
            allocate(mpbasis%radbasfn_mt(atoms%jmtd, nn + 1, 0:maxval(hybrid%lcutm1), atoms%ntype))
            mpbasis%radbasfn_mt(:,1:nn, :,:) = basmhlp(:,1:nn, :,:)
            deallocate(basmhlp)
         END IF

         ! add l = 0 function
         mpbasis%radbasfn_mt(:n_grid_pt, nn + 1, 0, itype) &
           = atoms%rmsh(:n_grid_pt, itype) / SQRT(atoms%rmsh(n_grid_pt, itype)**3/3)

         ! perform gram-schmidt orthonormalization
         DO i = nn, 1, -1
            DO j = i + 1, nn + 1
               mpbasis%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                  = mpbasis%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                     - intgrf( &
                                mpbasis%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                              * mpbasis%radbasfn_mt(:n_grid_pt, j, 0, itype), &
                              atoms, itype, gridf) &
                       * mpbasis%radbasfn_mt(:n_grid_pt, j, 0, itype)

            END DO

            ! renormalize
            norm = SQRT(intgrf(mpbasis%radbasfn_mt(:n_grid_pt, i, 0, itype)**2, atoms, itype, gridf))
            mpbasis%radbasfn_mt(:n_grid_pt, i, 0, itype) &
               = mpbasis%radbasfn_mt(:n_grid_pt, i, 0, itype) / norm
         END DO
         nn = nn + 1
         mpbasis%num_radbasfn(l, itype) = nn
      END IF
   end subroutine mpbasis_add_l0_fun

   subroutine mpbasis_reduce_linear_dep(mpbasis, atoms, mpi, hybrid, l, itype, gridf)
      use m_types_setup
      use m_types_hybrid
      use m_types_mpi
      implicit none
      class(t_mpbasis)              :: mpbasis
      type(t_atoms), intent(in)     :: atoms
      type(t_mpi), intent(in)       :: mpi
      type(t_hybrid), intent(in)    :: hybrid
      integer, intent(in)           :: l, itype

      REAL, ALLOCATABLE             :: olap(:,:), eig(:), eigv(:,:)
      REAL                          :: gridf(:,:)
      integer                       :: full_n_radbasfn, n_grid_pt

      full_n_radbasfn = mpbasis%num_radbasfn(l, itype)
      n_grid_pt = atoms%jri(itype)

      ! Calculate overlap
      call mpbasis%calc_olap_radbasfn(atoms, l, itype, gridf, olap)

      ! Diagonalize
      call mpbasis_diagonialize_olap(olap, eig, eigv)

      call mpbasis%filter_radbasfn(l, itype, full_n_radbasfn, eig, eigv)

      call mpbasis%trafo_to_orthonorm_bas(full_n_radbasfn, n_grid_pt, l, itype, eig, eigv)

      ! Add constant function to l=0 basis and then do a Gram-Schmidt orthonormalization
      call mpbasis%add_l0_fun(atoms, hybrid, n_grid_pt, l, itype, gridf)

      ! Check orthonormality of product basis
      call mpbasis%check_orthonormality(atoms, mpi, l, itype, gridf)

      deallocate(olap, eigv, eig)
   end subroutine
end module m_types_mpbasis
