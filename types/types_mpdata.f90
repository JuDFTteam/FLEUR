module m_types_mpdata
   implicit none

   type t_mpdata
      integer, allocatable   :: g(:, :) ! (3, num_gpts)
      integer, allocatable   :: n_g(:) ! (ik)
      integer, allocatable   :: gptm_ptr(:, :) ! (ig, ik)
      real                   :: g_cutoff
      integer, allocatable   :: num_radbasfn(:, :) !(l,itype)
      real, allocatable      :: radbasfn_mt(:,:,:,:) !(jri,n,l,itype)
      real                   :: linear_dep_tol  !only read in
      INTEGER, ALLOCATABLE   :: num_radfun_per_l(:,:) !(l,itype)

      integer, allocatable   :: l1(:,:,:) !(n, l, itype)
      integer, allocatable   :: l2(:,:,:) !(n, l, itype)
      integer, allocatable   :: n1(:,:,:) !(n, l, itype)
      integer, allocatable   :: n2(:,:,:) !(n, l, itype)
   CONTAINS
      procedure :: num_gpts               => mpdata_num_gpts
      procedure :: gen_gvec               => mpdata_gen_gvec
      procedure :: check_orthonormality   => mpdata_check_orthonormality
      procedure :: check_radbasfn         => mpdata_check_radbasfn
      procedure :: calc_olap_radbasfn     => mpdata_calc_olap_radbasfn
      procedure :: filter_radbasfn        => mpdata_filter_radbasfn
      procedure :: trafo_to_orthonorm_bas => mpdata_trafo_to_orthonorm_bas
      procedure :: add_l0_fun             => mpdata_add_l0_fun
      procedure :: reduce_linear_dep      => mpdata_reduce_linear_dep
      procedure :: normalize              => mpdata_normalize
      procedure :: set_nl                 => mpdata_set_nl
      procedure :: free                   => mpdata_free
      procedure :: init                   => mpdata_init
      procedure :: set_num_radfun_per_l   => set_num_radfun_per_l_mpdata
      !generic   :: write(unformatted) => write_mpdata
   end type t_mpdata
contains
   function mpdata_num_gpts(mpdata)
      implicit NONE
      class(t_mpdata), intent(in) :: mpdata
      integer    :: mpdata_num_gpts

      mpdata_num_gpts = size(mpdata%g, dim=2)
   end function mpdata_num_gpts

   subroutine mpdata_gen_gvec(mpdata, cell, kpts, mpi)
      use m_types_setup
      use m_types_kpts
      use m_types_mpi
      use m_intgrf, only: intgrf_init, intgrf
      use m_rorder, only: rorderpf
      implicit NONE
      class(t_mpdata), intent(inout) :: mpdata
      type(t_cell), intent(in)       :: cell
      type(t_kpts), intent(in)       :: kpts
      type(t_mpi), intent(in)        :: mpi

      integer :: i, n, n1, n2, divconq
      integer :: x, y, z, ikpt, igpt
      integer :: g(3)
      real    :: longest_k, kvec(3)
      logical :: l_found_new_gpt, l_found_kg_in_sphere

      integer, allocatable            ::  unsrt_pgptm(:, :) ! unsorted pointers to g vectors
      real, allocatable               ::  length_kg(:, :) ! length of the vectors k + G
      integer, allocatable            ::  ptr(:)

      allocate(mpdata%n_g(kpts%nkptf))

      mpdata%n_g = 0
      i = 0
      n = -1

      longest_k = MAXVAL([(norm2(MATMUL(kpts%bkf(:, ikpt), cell%bmat)), ikpt=1, kpts%nkptf)])

      ! a first run for the determination of the dimensions of the fields g,pgptm

      do
         n = n + 1
         l_found_new_gpt = .FALSE.
         do x = -n, n
            n1 = n - ABS(x)
            do y = -n1, n1
               n2 = n1 - ABS(y)
               do z = -n2, n2, MAX(2*n2, 1)
                  g = [x, y, z]
                  if ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpdata%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  do ikpt = 1, kpts%nkptf
                     if (norm2(MATMUL(kpts%bkf(:, ikpt) + g, cell%bmat)) <= mpdata%g_cutoff) THEN
                        if (.NOT. l_found_kg_in_sphere) THEN
                           i = i + 1
                           l_found_kg_in_sphere = .TRUE.
                        END if

                        mpdata%n_g(ikpt) = mpdata%n_g(ikpt) + 1
                        l_found_new_gpt = .TRUE.
                     END if
                  enddo ! k-loop
               enddo
            enddo
         enddo
         if (.NOT. l_found_new_gpt) EXIT
      enddo

      allocate(mpdata%g(3, i)) ! i = gptmd
      allocate(mpdata%gptm_ptr(maxval(mpdata%n_g), kpts%nkptf))

      ! allocate and initialize arrays needed for G vector ordering
      allocate(unsrt_pgptm(maxval(mpdata%n_g), kpts%nkptf))
      allocate(length_kG(maxval(mpdata%n_g), kpts%nkptf))

      mpdata%g = 0
      mpdata%gptm_ptr = 0
      mpdata%n_g = 0

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
                  if ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpdata%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  do ikpt = 1, kpts%nkptf
                     kvec = kpts%bkf(:, ikpt)

                     if (norm2(MATMUL(kvec + g, cell%bmat)) <= mpdata%g_cutoff) THEN
                        if (.NOT. l_found_kg_in_sphere) THEN
                           i = i + 1
                           mpdata%g(:, i) = g
                           l_found_kg_in_sphere = .TRUE.
                        END if

                        mpdata%n_g(ikpt) = mpdata%n_g(ikpt) + 1
                        l_found_new_gpt = .TRUE.

                        ! Position of the vector is saved as pointer
                        unsrt_pgptm(mpdata%n_g(ikpt), ikpt) = i
                        ! Save length of vector k + G for array sorting
                        length_kG(mpdata%n_g(ikpt), ikpt) = norm2(MATMUL(kvec + g, cell%bmat))
                     END if
                  enddo
               enddo
            enddo
         enddo
         if (.NOT. l_found_new_gpt) EXIT
      enddo

      ! Sort pointers in array, so that shortest |k+G| comes first
      do ikpt = 1, kpts%nkptf
         allocate(ptr(mpdata%n_g(ikpt)))
         ! Divide and conquer algorithm for arrays > 1000 entries
         divconq = MAX(0, INT(1.443*LOG(0.001*mpdata%n_g(ikpt))))
         ! create pointers which correspond to a sorted array
         CALL rorderpf(ptr, length_kG(1:mpdata%n_g(ikpt), ikpt), mpdata%n_g(ikpt), divconq)
         ! rearrange old pointers
         do igpt = 1, mpdata%n_g(ikpt)
            mpdata%gptm_ptr(igpt, ikpt) = unsrt_pgptm(ptr(igpt), ikpt)
         enddo
         deallocate(ptr)
      enddo

      if (mpi%irank == 0) THEN
         WRITE (6, '(/A)') 'Mixed basis'
         WRITE (6, '(A,I5)') 'Number of unique G-vectors: ', mpdata%num_gpts()
         WRITE (6, *)
         WRITE (6, '(3x,A)') 'IR Plane-wave basis with cutoff of gcutm (mpdata%g_cutoff/2*input%rkmax):'
         WRITE (6, '(5x,A,I5)') 'Maximal number of G-vectors:', maxval(mpdata%n_g)
      END if
   end subroutine mpdata_gen_gvec

   subroutine mpdata_check_orthonormality(mpdata, atoms, mpi, l, itype, gridf)
      use m_intgrf, only: intgrf
      use m_types_setup
      use m_types_mpi
      use m_judft
      implicit none
      class(t_mpdata)          :: mpdata
      type(t_atoms), intent(in) :: atoms
      type(t_mpi), intent(in)   :: mpi
      integer, intent(in)       :: itype, l
      real, intent(in)          :: gridf(:, :)

      real, allocatable :: olap(:, :)
      integer           :: i, n_radbasfn, err_loc(2)

      CHARACTER, PARAMETER :: lchar(0:38) = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', &
                                             'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', &
                                             'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x']

      call timestart("check mpdata orthonormality")

      ! calculate overlap matrix
      call mpdata%calc_olap_radbasfn(atoms, l, itype, gridf, olap)

      !subtract identity-matrix
      do i = 1, size(olap, dim=1)
         olap(i, i) = olap(i, i) - 1.0
      enddo

      ! check if (olap - identity) is zero-matrix
      if (norm2(olap) > 1e-6) then
         if (mpi%irank == 0) THEN
            err_loc = maxloc(abs(olap))
            WRITE (6, '(A)') 'mixedbasis: Bad orthonormality of ' &
               //lchar(l)//'-mpdatauct basis. Increase tolerance.'
            WRITE (6, '(12X,A,F9.6,A,2(I3.3,A))') 'Deviation of', &
               maxval(abs(olap)), ' occurred for (', &
               err_loc(1), ',', err_loc(2), ')-overlap.'
         endif
         CALL judft_error("Bad orthonormality of mpdata", &
                          hint='Increase tolerance', &
                          calledby='mixedbasis%check_orthonormality')
      endif

      if (mpi%irank == 0) THEN
         n_radbasfn = mpdata%num_radbasfn(l, itype)
         WRITE (6, '(6X,A,I4,''   ('',ES8.1,'' )'')') &
            lchar(l)//':', n_radbasfn, norm2(olap)/n_radbasfn
      END if
      call timestop("check mpdata orthonormality")
   end subroutine mpdata_check_orthonormality

   subroutine mpdata_check_radbasfn(mpdata, atoms, hybinp)
      use m_judft
      use m_types_hybinp
      use m_types_setup
      implicit none
      class(t_mpdata), intent(in) :: mpdata
      type(t_atoms), intent(in)    :: atoms
      type(t_hybinp), intent(in)   :: hybinp

      integer :: itype

      do itype = 1, atoms%ntype
         if (ANY(mpdata%num_radbasfn(0:hybinp%lcutm1(itype), itype) == 0)) THEN
            call judft_error('any mpdata%num_radbasfn eq 0', calledby='mixedbasis')
         endif
      enddo
   end subroutine mpdata_check_radbasfn

   subroutine mpdata_calc_olap_radbasfn(mpdata, atoms, l, itype, gridf, olap)
      use m_intgrf, only: intgrf
      use m_types_setup
      use m_judft
      implicit NONE
      class(t_mpdata), intent(in)       :: mpdata
      type(t_atoms), intent(in)          :: atoms
      integer, intent(in)                :: l, itype
      real, intent(in)                   :: gridf(:, :)
      real, intent(inout), allocatable   :: olap(:, :)

      integer  :: n1, n2, n_radbasfn

      call timestart("calc mpdata overlap")

      n_radbasfn = mpdata%num_radbasfn(l, itype)
      if (allocated(olap)) then
         if (any(shape(olap) /= n_radbasfn)) then
            deallocate(olap)
         endif
      endif
      if (.not. allocated(olap)) allocate(olap(n_radbasfn, n_radbasfn), source=0.0)

      do n2 = 1, n_radbasfn
         do n1 = 1, n2
            olap(n1, n2) = intgrf(mpdata%radbasfn_mt(:, n1, l, itype)*mpdata%radbasfn_mt(:, n2, l, itype), &
                                  atoms, itype, gridf)
            if(isnan(olap(n1,n2))) then
               write (*,*) "nan at", n1, n2
            endif
            olap(n2, n1) = olap(n1, n2)
         END do
      END do
      if(any(isnan(olap))) call juDFT_error("Mixed-product basis olap is nan")
      call timestop("calc mpdata overlap")
   end subroutine mpdata_calc_olap_radbasfn

   subroutine mpdata_filter_radbasfn(mpdata, l, itype, n_radbasfn, eig, eigv)
      ! Get rid of linear dependencies (eigenvalue <= mpdata%linear_dep_tol)
      use m_judft
      implicit none
      class(t_mpdata), intent(inout)       :: mpdata
      integer, intent(in)                   :: l, itype, n_radbasfn
      real, intent(inout)                   :: eig(:), eigv(:, :)

      integer              :: num_radbasfn, i_bas
      integer, allocatable :: remaining_basfn(:)

      call timestart("filer mpdata")
      allocate(remaining_basfn(n_radbasfn), source=1)
      num_radbasfn = 0

      do i_bas = 1, mpdata%num_radbasfn(l, itype)
         if (eig(i_bas) > mpdata%linear_dep_tol) THEN
            num_radbasfn = num_radbasfn + 1
            remaining_basfn(num_radbasfn) = i_bas
         END if
      END do

      mpdata%num_radbasfn(l, itype) = num_radbasfn
      eig = eig(remaining_basfn)
      eigv(:, :) = eigv(:, remaining_basfn)
      call timestop("filer mpdata")
   end subroutine mpdata_filter_radbasfn

   subroutine mpdata_diagonialize_olap(olap, eig_val, eig_vec)
      use m_judft
      implicit NONE
      real, intent(in)  :: olap(:, :)
      real, allocatable :: eig_val(:), eig_vec(:, :)

      integer              :: n, size_iwork, info
      real                 :: size_work
      integer, allocatable :: iwork(:)
      real, allocatable    :: work(:)

      call timestart("diagonalize overlap")
      if (size(olap, dim=1) /= size(olap, dim=2)) then
         call juDFT_error("only square matrices can be diagonalized")
      endif

      n = size(olap, dim=1)

      if (allocated(eig_val)) then
         if (size(eig_val) /= n) deallocate(eig_val)
      endif
      if (.not. allocated(eig_val)) allocate(eig_val(n))

      eig_vec = olap
      ! get sizes of work arrays
      call dsyevd('V', 'U', n, eig_vec, n, eig_val, &
                  size_work, -1, size_iwork, -1, info)
      if (info /= 0) call juDFT_error("diagonalization for size failed")

      allocate(work(int(size_work)))
      allocate(iwork(size_iwork))

      call dsyevd('V', 'U', n, eig_vec, n, eig_val, &
                  work, int(size_work), iwork, size_iwork, info)
      if (info /= 0) call juDFT_error("diagonalization failed")
      call timestop("diagonalize overlap")
   end subroutine mpdata_diagonialize_olap

   subroutine mpdata_trafo_to_orthonorm_bas(mpdata, full_n_radbasfn, n_grid_pt, l, itype, eig, eigv)
      use m_judft
      implicit NONE
      class(t_mpdata), intent(inout)  :: mpdata
      integer, intent(in)              :: full_n_radbasfn, n_grid_pt, l, itype
      real, intent(in)                 :: eig(:), eigv(:, :)

      integer :: nn, i

      call timestart("transform to reduced mpdata")
      ! reduced number of basis functions
      nn = mpdata%num_radbasfn(l, itype)

      do i = 1, n_grid_pt
         mpdata%radbasfn_mt(i, 1:nn, l, itype) &
            = MATMUL(mpdata%radbasfn_mt(i, 1:full_n_radbasfn, l, itype), eigv(:, 1:nn))/SQRT(eig(:nn))
      END do
      call timestop("transform to reduced mpdata")
   end subroutine mpdata_trafo_to_orthonorm_bas

   subroutine mpdata_add_l0_fun(mpdata, atoms, hybinp, n_grid_pt, l, itype, gridf)
      use m_types_setup
      use m_types_hybinp
      use m_intgrf, only: intgrf
      use m_judft
      implicit none
      class(t_mpdata), intent(inout) :: mpdata
      type(t_atoms), intent(in)       :: atoms
      type(t_hybinp), intent(in)      :: hybinp
      integer, intent(in)             :: n_grid_pt, l, itype
      real, intent(in)                :: gridf(:, :)

      integer                         :: i, j, nn
      real, allocatable               :: basmhlp(:, :, :, :)
      real                            :: norm

      call timestart("add l0 to mpdata")
      nn = mpdata%num_radbasfn(l, itype)
      if (l == 0) THEN

         ! Check if radbasfn_mt must be reallocated
         if (nn + 1 > SIZE(mpdata%radbasfn_mt, 2)) THEN
            allocate(basmhlp(atoms%jmtd, nn + 1, 0:maxval(hybinp%lcutm1), atoms%ntype))
            basmhlp(:, 1:nn, :, :) = mpdata%radbasfn_mt
            deallocate(mpdata%radbasfn_mt)
            allocate(mpdata%radbasfn_mt(atoms%jmtd, nn + 1, 0:maxval(hybinp%lcutm1), atoms%ntype))
            mpdata%radbasfn_mt(:, 1:nn, :, :) = basmhlp(:, 1:nn, :, :)
            deallocate(basmhlp)
         END if

         ! add l = 0 function
         mpdata%radbasfn_mt(:n_grid_pt, nn + 1, 0, itype) &
            = atoms%rmsh(:n_grid_pt, itype)/SQRT(atoms%rmsh(n_grid_pt, itype)**3/3)

         ! perform gram-schmidt orthonormalization
         do i = nn, 1, -1
            do j = i + 1, nn + 1
               mpdata%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                  = mpdata%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                    - intgrf( &
                    mpdata%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                    *mpdata%radbasfn_mt(:n_grid_pt, j, 0, itype), &
                    atoms, itype, gridf) &
                    *mpdata%radbasfn_mt(:n_grid_pt, j, 0, itype)

            END do

            ! renormalize
            norm = SQRT(intgrf(mpdata%radbasfn_mt(:n_grid_pt, i, 0, itype)**2, atoms, itype, gridf))
            mpdata%radbasfn_mt(:n_grid_pt, i, 0, itype) &
               = mpdata%radbasfn_mt(:n_grid_pt, i, 0, itype)/norm
         END do
         nn = nn + 1
         mpdata%num_radbasfn(l, itype) = nn
      END if
      call timestop("add l0 to mpdata")
   end subroutine mpdata_add_l0_fun

   subroutine mpdata_reduce_linear_dep(mpdata, atoms, mpi, hybinp, l, itype, gridf, iterHF)
      use m_types_setup
      use m_types_hybinp
      use m_types_mpi
      use m_judft
      implicit none
      class(t_mpdata)              :: mpdata
      type(t_atoms), intent(in)     :: atoms
      type(t_mpi), intent(in)       :: mpi
      type(t_hybinp), intent(in)    :: hybinp
      integer, intent(in)           :: l, itype, iterHF

      real, allocatable             :: olap(:, :), eig(:), eigv(:, :)
      real                          :: gridf(:, :)
      integer                       :: full_n_radbasfn, n_grid_pt

      call timestart("reduce lin. dep. mpdata")
      full_n_radbasfn = mpdata%num_radbasfn(l, itype)
      n_grid_pt = atoms%jri(itype)

      ! Calculate overlap
      call mpdata%calc_olap_radbasfn(atoms, l, itype, gridf, olap)

      ! Diagonalize
      call mpdata_diagonialize_olap(olap, eig, eigv)

      call mpdata%filter_radbasfn(l, itype, full_n_radbasfn, eig, eigv)

      call mpdata%trafo_to_orthonorm_bas(full_n_radbasfn, n_grid_pt, l, itype, eig, eigv)

      ! Add constant function to l=0 basis and then do a Gram-Schmidt orthonormalization
      call mpdata%add_l0_fun(atoms, hybinp, n_grid_pt, l, itype, gridf)

      ! Check orthonormality of mpdatauct basis
      call mpdata%check_orthonormality(atoms, mpi, l, itype, gridf)

      deallocate(olap, eigv, eig)
      call timestop("reduce lin. dep. mpdata")
   end subroutine

   subroutine mpdata_normalize(mpdata, atoms, hybinp, gridf)
      use m_intgrf, only: intgrf
      use m_types_hybinp
      use m_types_setup
      implicit NONE

      class(t_mpdata), intent(inout):: mpdata
      type(t_atoms), intent(in)      :: atoms
      type(t_hybinp), intent(in)     :: hybinp
      real, intent(in)               :: gridf(:, :)

      integer                        :: l, i_basfn, itype
      real                           :: norm

      do itype = 1, atoms%ntype
         do l = 0, hybinp%lcutm1(itype)
            do i_basfn = 1, mpdata%num_radbasfn(l, itype)
               norm = SQRT( &
                      intgrf(mpdata%radbasfn_mt(:, i_basfn, l, itype)**2, &
                             atoms, itype, gridf) &
                      )

               mpdata%radbasfn_mt(:atoms%jri(itype), i_basfn, l, itype) &
                  = mpdata%radbasfn_mt(:atoms%jri(itype), i_basfn, l, itype)/norm
            end do
         end do
      end do
   end subroutine mpdata_normalize

   subroutine mpdata_init(mpdata, hybinp, hybdat, atoms)
      use m_types_setup
      use m_types_hybinp
      use m_types_hybdat
      use m_judft
      implicit none
      class(t_mpdata)           :: mpdata
      type(t_hybinp), intent(in) :: hybinp
      type(t_hybdat), intent(in) :: hybdat
      type(t_atoms), intent(in)  :: atoms

      integer                    :: ok

      if(.not. allocated(mpdata%l1)) then
         allocate(mpdata%l1(hybdat%max_indx_p_1, 0:maxval(hybinp%lcutm1), atoms%ntype), stat=ok)
         if (ok /= 0) call judft_error('mpdata_init: failure allocation mpdata%l1')

         allocate(mpdata%l2, mold=mpdata%l1, stat=ok)
         if (ok /= 0) call judft_error('mpdata_init: failure allocation mpdata%l2')

         allocate(mpdata%n1, mold=mpdata%l1, stat=ok)
         if (ok /= 0) call judft_error('mpdata_init: failure allocation mpdata%n1')

         allocate(mpdata%n2, mold=mpdata%l1, stat=ok)
         if (ok /= 0) call judft_error('mpdata_init: failure allocation mpdata%n2')
      endif
   end subroutine mpdata_init

   subroutine mpdata_free(mpdata)
      implicit NONE
      class(t_mpdata)          :: mpdata

      if (allocated(mpdata%l1)) deallocate(mpdata%l1)
      if (allocated(mpdata%l2)) deallocate(mpdata%l2)
      if (allocated(mpdata%n1)) deallocate(mpdata%n1)
      if (allocated(mpdata%n2)) deallocate(mpdata%n2)
   end subroutine mpdata_free

   subroutine mpdata_set_nl(mpdata, n, l, itype, n1, l1, n2, l2)
      implicit NONE
      class(t_mpdata)    :: mpdata
      integer, intent(in)  :: n, l, itype
      integer, intent(out) :: n1, l1, n2, l2

      l1 = mpdata%l1(n, l, itype) !
      l2 = mpdata%l2(n, l, itype) ! current basis-function mpdatauct
      n1 = mpdata%n1(n, l, itype) ! = bas(:,n1,l1,itype)*bas(:,n2,l2,itype) = b1*b2
      n2 = mpdata%n2(n, l, itype) !

   end subroutine mpdata_set_nl

   subroutine set_num_radfun_per_l_mpdata(mpdata, atoms)
      use m_types_setup
      implicit NONE
      class(t_mpdata) :: mpdata
      type(t_atoms)   :: atoms
      integer :: itype, ilo

      ! there is always at least two: u and u_dot
      mpdata%num_radfun_per_l = 2
      DO itype = 1, atoms%ntype
         DO ilo = 1, atoms%nlo(itype)
            mpdata%num_radfun_per_l(atoms%llo(ilo, itype), itype) &
              = mpdata%num_radfun_per_l(atoms%llo(ilo, itype), itype) + 1
         END DO
      END DO
   end subroutine set_num_radfun_per_l_mpdata
end module m_types_mpdata
