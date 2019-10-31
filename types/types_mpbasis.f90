module m_types_mpbasis
   implicit none

   type t_mpbasis
      INTEGER, ALLOCATABLE   ::  gptm(:,:) ! (3, num_gpts)
      INTEGER, ALLOCATABLE   ::  ngptm(:)
      INTEGER, ALLOCATABLE   ::  gptm_ptr(:,:)
      REAL                   ::  g_cutoff
   CONTAINS
      procedure :: num_gpts => mpbasis_num_gpts
      procedure :: gen_gvec => mpbasis_gen_gvec
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

      DO
         n = n + 1
         l_found_new_gpt = .FALSE.
         DO x = -n, n
            n1 = n - ABS(x)
            DO y = -n1, n1
               n2 = n1 - ABS(y)
               DO z = -n2, n2, MAX(2*n2, 1)
                  g = [x, y, z]
                  IF ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpbasis%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  DO ikpt = 1, kpts%nkptf
                     IF (norm2(MATMUL(kpts%bkf(:,ikpt) + g, cell%bmat)) <= mpbasis%g_cutoff) THEN
                        IF (.NOT. l_found_kg_in_sphere) THEN
                           i = i + 1
                           l_found_kg_in_sphere = .TRUE.
                        END IF

                        mpbasis%ngptm(ikpt) = mpbasis%ngptm(ikpt) + 1
                        l_found_new_gpt = .TRUE.
                     END IF
                  END DO ! k-loop
               END DO
            END DO
         END DO
         IF (.NOT. l_found_new_gpt) EXIT
      END DO

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

      DO
         n = n + 1
         l_found_new_gpt = .FALSE.
         DO x = -n, n
            n1 = n - ABS(x)
            DO y = -n1, n1
               n2 = n1 - ABS(y)
               DO z = -n2, n2, MAX(2*n2, 1)
                  g = [x, y, z]
                  IF ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpbasis%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  DO ikpt = 1, kpts%nkptf
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
                  END DO
               END DO
            END DO
         END DO
         IF (.NOT. l_found_new_gpt) EXIT
      END DO

      ! Sort pointers in array, so that shortest |k+G| comes first
      DO ikpt = 1, kpts%nkptf
         allocate(ptr(mpbasis%ngptm(ikpt)))
         ! Divide and conquer algorithm for arrays > 1000 entries
         divconq = MAX(0, INT(1.443*LOG(0.001*mpbasis%ngptm(ikpt))))
         ! create pointers which correspond to a sorted array
         CALL rorderpf(ptr, length_kG(1:mpbasis%ngptm(ikpt), ikpt), mpbasis%ngptm(ikpt), divconq)
         ! rearrange old pointers
         DO igpt = 1, mpbasis%ngptm(ikpt)
            mpbasis%gptm_ptr(igpt, ikpt) = unsrt_pgptm(ptr(igpt), ikpt)
         END DO
         deallocate(ptr)
      END DO

      IF (mpi%irank == 0) THEN
         WRITE (6, '(/A)') 'Mixed basis'
         WRITE (6, '(A,I5)') 'Number of unique G-vectors: ', mpbasis%num_gpts()
         WRITE (6, *)
         WRITE (6, '(3x,A)') 'IR Plane-wave basis with cutoff of gcutm (mpbasis%g_cutoff/2*input%rkmax):'
         WRITE (6, '(5x,A,I5)') 'Maximal number of G-vectors:', maxval(mpbasis%ngptm)
      END IF
   end subroutine mpbasis_gen_gvec
end module m_types_mpbasis
