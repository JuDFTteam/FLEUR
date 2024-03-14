MODULE m_types_hybdat
   use m_types_usdus
   use m_types_mat
   use m_types_coul
   use m_types_eigvec
   use m_io_matrix
   use m_judft
   use m_types_misc
   use m_types_fleurinput
   use m_types_mpi
   use m_constants
      
#ifdef CPP_MPI
   use mpi
#endif
   IMPLICIT NONE
   private 
   TYPE,public:: t_hybdat
      COMPLEX, ALLOCATABLE   :: stepfunc(:, :, :)
      INTEGER, ALLOCATABLE   :: lmaxc(:)
      INTEGER, ALLOCATABLE   :: nbands(:,:) ! nkptf, jsp
      INTEGER, ALLOCATABLE   :: nbasm(:)
      INTEGER, ALLOCATABLE   :: nindxc(:, :)
      INTEGER, ALLOCATABLE   :: nindxp1(:, :)
      INTEGER, ALLOCATABLE   :: nobd(:, :)
      INTEGER, ALLOCATABLE   :: pntgpt(:, :, :, :)
      INTEGER, ALLOCATABLE   :: pntgptd(:)
      INTEGER                :: eig_id = -1
      INTEGER                :: lmaxcd, maxindxc
      INTEGER                :: maxfac
      INTEGER                :: maxlmindx = -1
      INTEGER                :: nbasp = -1
      integer                :: max_q = -1
      LOGICAL                :: l_addhf = .false.
      LOGICAL                :: l_calhf = .false.
      LOGICAL                :: l_print_iob_splitting = .True.
      LOGICAL                :: l_subvxc = .false.
      REAL, ALLOCATABLE      :: bas1(:, :, :, :), bas2(:, :, :, :)
      REAL, ALLOCATABLE      :: bas1_MT(:, :, :), drbas1_MT(:, :, :)
      REAL, ALLOCATABLE      :: core1(:, :, :, :), core2(:, :, :, :)
      REAL, ALLOCATABLE      :: div_vv(:, :, :)
      REAL, ALLOCATABLE      :: eig_c(:, :, :)
      REAL, ALLOCATABLE      :: gauntarr(:, :, :, :, :, :)
      REAL, ALLOCATABLE      :: gridf(:, :)
      REAL, ALLOCATABLE      :: prodm(:, :, :, :)
      REAL, ALLOCATABLE      :: sfac(:), fac(:)
      ! coulomb matrix stuff
      type(t_coul), allocatable   :: coul(:)

      type(t_usdus)               :: usdus
      class(t_mat), allocatable   :: v_x(:, :) ! nkpt, jsp
      type(t_eigvec), allocatable :: zmat(:, :) ! nkpt, jsp
      type(t_results)             :: results
   contains
      procedure :: set_stepfunction => set_stepfunction
      procedure :: free             => free_hybdat
      procedure :: allocate         => allocate_hybdat
      procedure :: set_nobd         => set_nobd_hybdat
      procedure :: set_nbands       => set_nbands_hybdat
      procedure :: set_maxlmindx    => set_maxlmindx_hybdat
   END TYPE t_hybdat
   public:: gptnorm
contains
   subroutine set_maxlmindx_hybdat(hybdat, atoms, num_radfun_per_l)
      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_atoms), intent(in)      :: atoms
      integer, intent(in)            :: num_radfun_per_l(0:,:)

      integer :: itype, l, summe

      hybdat%maxlmindx = 0
      do itype = 1,atoms%ntype
         summe = 0 
         do l = 0, atoms%lmax(itype)
            summe = summe + (2*l + 1) * num_radfun_per_l(l, itype)
         enddo
         hybdat%maxlmindx = max(hybdat%maxlmindx, summe)
      enddo
   end subroutine set_maxlmindx_hybdat

   subroutine set_nobd_hybdat(hybdat, fi, results)
      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in) :: fi
      type(t_results), intent(in)    :: results

      integer :: jsp, nk, i

      if (.not. allocated(hybdat%nobd)) then
         allocate (hybdat%nobd(fi%kpts%nkptf, fi%input%jspins))
      endif
      hybdat%nobd = 0

      do jsp = 1, fi%input%jspins
         DO nk = 1, fi%kpts%nkpt
            DO i = 1, results%neig(nk, jsp)
               IF (results%w_iks(i, nk, jsp) > 0.0) then
                  hybdat%nobd(nk, jsp) = hybdat%nobd(nk, jsp) + 1
               endif
            END DO
         enddo
      enddo

      do jsp = 1, fi%input%jspins
         DO nk = 1, fi%kpts%nkptf
            i = fi%kpts%bkp(nk)
            hybdat%nobd(nk, jsp) = hybdat%nobd(i, jsp)
         END DO
      enddo
   end subroutine set_nobd_hybdat

   subroutine set_nbands_hybdat(hybdat, fi, fmpi, results)
      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in) :: fi
      type(t_mpi), intent(in)        :: fmpi
      type(t_results), intent(in)    :: results


      integer :: i, j, nk, jsp

      INTEGER :: degenerat(merge(fi%input%neig*2,fi%input%neig,fi%noco%l_soc) + 1, fi%kpts%nkpt)
      !determine degenerate states at each k-point
      !
      ! degenerat(i) =1  band i  is not degenerat ,
      ! degenerat(i) =j  band i  has j-1 degenart states ( i, i+1, ..., i+j)
      ! degenerat(i) =0  band i  is  degenerat, but is not the lowest band
      !                  of the group of degenerate states

      call timestart("degenerate treatment")
      IF (fmpi%irank == 0) THEN
         WRITE (oUnit, *)
         WRITE (oUnit, '(A)') "   k-point      |   number of occupied bands  |   maximal number of bands"
      END IF

      do jsp = 1,fi%input%jspins
         degenerat = 1
         DO nk = 1, fi%kpts%nkpt
            DO i = 1, results%neig(nk, jsp)
               DO j = i + 1, results%neig(nk, jsp)
                  IF (ABS(results%eig(i, nk, jsp) - results%eig(j, nk, jsp)) < 1E-07) THEN !0.015
                     degenerat(i, nk) = degenerat(i, nk) + 1
                  END IF
               END DO
            END DO

            DO i = 1, results%neig(nk, jsp)
               IF ((degenerat(i, nk) /= 1) .OR. (degenerat(i, nk) /= 0)) degenerat(i + 1:i + degenerat(i, nk) - 1, nk) = 0
            END DO

            ! set the size of the exchange matrix in the space of the wavefunctions

            hybdat%nbands(nk,jsp) = fi%hybinp%bands1
            IF (hybdat%nbands(nk,jsp) > results%neig(nk, jsp)) THEN
               IF (fmpi%irank == 0) THEN
                  WRITE (*, *) ' maximum for hybdat%nbands is', results%neig(nk, jsp)
                  WRITE (*, *) ' increase energy window to obtain enough eigenvalues'
                  WRITE (*, *) ' set hybdat%nbands equal to results%neig'
               END IF
               hybdat%nbands(nk,jsp) = results%neig(nk, jsp)
            END IF

            DO i = hybdat%nbands(nk,jsp) - 1, 1, -1
               IF ((degenerat(i, nk) >= 1) .AND. (degenerat(i, nk) + i - 1 /= hybdat%nbands(nk,jsp))) THEN
                  hybdat%nbands(nk,jsp) = i + degenerat(i, nk) - 1
                  EXIT
               END IF
            END DO

            IF (hybdat%nobd(nk, jsp) > hybdat%nbands(nk,jsp)) THEN
               WRITE (*, *) 'k-point: ', nk
               WRITE (*, *) 'number of bands:          ', hybdat%nbands(nk,jsp)
               WRITE (*, *) 'number of occupied bands: ', hybdat%nobd(nk, jsp)
               CALL judft_warn("More occupied bands than total no of bands!?")
               hybdat%nbands(nk,jsp) = hybdat%nobd(nk, jsp)
            END IF
         END DO
         
         ! spread nbands from IBZ to whole BZ
         DO nk = fi%kpts%nkpt+1, fi%kpts%nkptf
            i = fi%kpts%bkp(nk)
            hybdat%nbands(nk,jsp) = hybdat%nbands(i,jsp)
         END DO
      enddo
      call timestop("degenerate treatment")
   end subroutine set_nbands_hybdat

   subroutine allocate_hybdat(hybdat, fi, num_radfun_per_l)

      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: num_radfun_per_l(:, :)
      integer                        :: ok(12)

      ok = -1
      allocate (hybdat%lmaxc(fi%atoms%ntype), &
                stat=ok(1), source=0)
      allocate (hybdat%bas1(fi%atoms%jmtd, maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
                stat=ok(2), source=0.0)
      allocate (hybdat%bas2(fi%atoms%jmtd, maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
                stat=ok(3), source=0.0)
      allocate (hybdat%bas1_MT(maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
                stat=ok(4), source=0.0)
      allocate (hybdat%drbas1_MT(maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
                stat=ok(5), source=0.0)

      ! core allocs
      allocate (hybdat%nindxc(0:hybdat%lmaxcd, fi%atoms%ntype), &
                stat=ok(6), source=0)
      allocate (hybdat%core1(fi%atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
                stat=ok(7), source=0.0)
      allocate (hybdat%core2(fi%atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
                stat=ok(8), source=0.0)
      allocate (hybdat%eig_c(hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
                stat=ok(9), source=0.0)

      allocate (hybdat%fac(0:hybdat%maxfac), stat=ok(10), source=0.0)
      allocate (hybdat%sfac(0:hybdat%maxfac), stat=ok(11), source=0.0)

      ALLOCATE (hybdat%gauntarr(2, 0:fi%atoms%lmaxd, 0:fi%atoms%lmaxd, 0:maxval(fi%hybinp%lcutm1), &
                                -fi%atoms%lmaxd:fi%atoms%lmaxd, -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1)), &
                stat=ok(12), source=0.0)

      if (any(ok /= 0)) then
         write (*, *) "allocation of hybdat failed. Error in array no.:"
         write (*, *) maxloc(abs(ok))
         call juDFT_error("allocation of hybdat failed. Error in array no is the outfile")
      endif
   end subroutine allocate_hybdat

   subroutine free_hybdat(hybdat)
      implicit none
      class(t_hybdat), intent(inout) :: hybdat

      if (allocated(hybdat%lmaxc)) deallocate (hybdat%lmaxc)
      if (allocated(hybdat%bas1)) deallocate (hybdat%bas1)
      if (allocated(hybdat%bas2)) deallocate (hybdat%bas2)
      if (allocated(hybdat%bas1_MT)) deallocate (hybdat%bas1_MT)
      if (allocated(hybdat%drbas1_MT)) deallocate (hybdat%drbas1_MT)
      if (allocated(hybdat%nindxc)) deallocate (hybdat%nindxc)
      if (allocated(hybdat%core1)) deallocate (hybdat%core1)
      if (allocated(hybdat%core2)) deallocate (hybdat%core2)
      if (allocated(hybdat%eig_c)) deallocate (hybdat%eig_c)
      if (allocated(hybdat%fac)) deallocate (hybdat%fac)
      if (allocated(hybdat%sfac)) deallocate (hybdat%sfac)
      if (allocated(hybdat%gauntarr)) deallocate (hybdat%gauntarr)
   end subroutine free_hybdat

   subroutine set_stepfunction(hybdat, cell, atoms, g, svol)
  
      implicit none
      class(t_hybdat), INTENT(INOUT) :: hybdat
      type(t_cell), INTENT(in)    :: cell
      type(t_atoms), INTENT(in)    :: atoms
      integer, INTENT(in)    :: g(3)
      real, INTENT(in)    :: svol
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
                  hybdat%stepfunc(i, j, k) = stepfunction(cell, atoms, [i, j, k])/svol
               END DO
            END DO
         END DO
         call timestop("setup stepfunction")
      endif

   end subroutine set_stepfunction

   !private subroutine
   FUNCTION stepfunction(cell, atoms, g)
   
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

   subroutine calc_matrix_slots(l_real, mtx_dim, slots_per_mtx, col_in_slot)
      implicit none
      logical, intent(in)  :: l_real
      integer, intent(in)  :: mtx_dim
      integer, intent(out) :: slots_per_mtx, col_in_slot

      integer(8)            :: mtx_size, type_size
      integer(8), parameter :: max_bytes = huge(slots_per_mtx) - 1
      integer  :: i
      type_size = merge(8, 16, l_real)

      ! avoid int32 overflow
      mtx_size = type_size*mtx_dim
      mtx_size = mtx_size*mtx_dim

      i = 1
      do while ((mtx_size/i >= max_bytes) .or. mod(mtx_dim, i) /= 0)
         i = i + 1
      enddo

      slots_per_mtx = i
      col_in_slot = mtx_dim/i
   end subroutine calc_matrix_slots
END MODULE m_types_hybdat
