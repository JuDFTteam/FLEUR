module m_wavefproducts_aux
   use m_types_fftGrid
CONTAINS
   subroutine wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                   ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)
      !$ use omp_lib
      use m_types
      use m_constants
      use m_judft
      use m_fft_interface
      use m_io_hybrid
      use m_juDFT
#ifdef CPP_MPI
      use mpi
#endif
      implicit NONE
      type(t_fleurinput), intent(in)  :: fi
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_mpdata), intent(in)      :: mpdata
      TYPE(t_hybdat), INTENT(IN)      :: hybdat
      type(t_stars), intent(in)       :: stars
      type(t_mat), intent(in)         :: z_k
      type(t_mat), intent(inout)      :: z_kqpt_p, cprod
      !     - scalars -
      INTEGER, INTENT(IN)      ::  ik, iq, jsp, g_t(3), bandoi, bandof
      INTEGER, INTENT(IN)      ::  ikqpt
      !     - arrays -
      complex, intent(inout)    :: c_phase_kqpt(hybdat%nbands(ikqpt,jsp))

      complex, allocatable  :: prod(:), psi_k(:, :)

      type(t_mat)     :: z_kqpt
      type(t_lapw)    :: lapw_ikqpt
      type(t_mat)     :: psi_kqpt
      type(t_fft)     :: fft
      type(t_fftgrid) :: stepf

      integer :: g(3), igptm, iob, n_omp
      integer :: ok, nbasfcn, psize, iband, ierr, i
      integer, allocatable :: band_list(:)
      real    :: inv_vol, gcutoff
      complex :: inv_gridlen

      logical :: real_warned

      real_warned = .False.

      call timestart("wavef_IS_FFT")
      gcutoff = (2*fi%input%rkmax + fi%mpinp%g_cutoff) * fi%hybinp%fftcut
      inv_vol = 1/sqrt(fi%cell%omtil)
      psize = bandof - bandoi + 1
      !this is for the exact result. Christoph recommend 2*gmax+gcutm for later
      if (2*fi%input%rkmax + fi%mpinp%g_cutoff > fi%input%gmax) then
         write (*, *) "WARNING: not accurate enough: 2*kmax+gcutm >= fi%input%gmax"
         !call juDFT_error("not accurate enough: 2*kmax+gcutm >= fi%input%gmax")
      endif

      call stepf%init(fi%cell, fi%sym, gcutoff)
      call stepf%putFieldOnGrid(stars, stars%ustep)
      stepf%grid = stepf%grid * inv_vol
      inv_gridlen = 1.0/stepf%gridLength

      call fft%init(stepf%dimensions, .false.)
      call fft%exec(stepf%grid)
      call fft%free()

      CALL lapw_ikqpt%init(fi, nococonv, ikqpt)

      nbasfcn = lapw_ikqpt%hyb_num_bas_fun(fi)
      call z_kqpt%alloc(z_k%l_real, nbasfcn, psize)
      call z_kqpt_p%init(z_kqpt)

      band_list = [(i, i=bandoi, bandof)]
      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, fi%input, ikqpt, jsp, z_kqpt, &
                  c_phase=c_phase_kqpt, parent_z=z_kqpt_p, list=band_list)

      call psi_kqpt%alloc(.false., stepf%gridLength, psize)

      call timestart("1st wavef2rs")
      call wavef2rs(fi, lapw_ikqpt, z_kqpt, gcutoff, 1, psize, jsp, psi_kqpt%data_c)
      call timestop("1st wavef2rs")

      call timestart("Big OMP loop")
      !$OMP PARALLEL default(none) &
      !$OMP private(iband, iob, g, igptm, prod, psi_k, ok, fft) &
      !$OMP shared(hybdat, psi_kqpt, cprod,  mpdata, iq, g_t, psize, gcutoff, inv_gridlen)&
      !$OMP shared(jsp, z_k, stars, lapw, fi, inv_vol, ik, real_warned, n_omp, bandoi, stepf)

      allocate (prod(0:stepf%gridLength - 1), stat=ok)
      if (ok /= 0) call juDFT_error("can't alloc prod")
      allocate (psi_k(0:stepf%gridLength - 1, 1), stat=ok)
      if (ok /= 0) call juDFT_error("can't alloc psi_k")

      call fft%init(stepf%dimensions, .true.)
      !$OMP DO
      do iband = 1, hybdat%nbands(ik,jsp)
         ! call timestart("loop wavef2rs")
         call wavef2rs(fi, lapw, z_k, gcutoff, iband, iband, jsp, psi_k)
         ! call timestop("loop wavef2rs")
         psi_k(:, 1) = conjg(psi_k(:, 1)) * stepf%grid

         do iob = 1, psize
            prod = psi_k(:, 1)*psi_kqpt%data_c(:, iob)
            call fft%exec(prod)
            if (cprod%l_real) then
               if (any(abs(aimag(prod)) > 1e-8) .and. (.not. real_warned)) then
                  write (*, *) "Imag part non-zero in is_fft maxval(abs(aimag(prod)))) = "// &
                     float2str(maxval(abs(aimag(prod))))
                  real_warned = .True.
               endif
            endif

            ! we still have to devide by the number of mesh points
            call zscal(stepf%gridLength, inv_gridlen, prod, 1)

            if (cprod%l_real) then
               DO igptm = 1, mpdata%n_g(iq)
                  g = mpdata%g(:, mpdata%gptm_ptr(igptm, iq)) - g_t
                  cprod%data_r(hybdat%nbasp + igptm, iob + (iband - 1)*psize) = real(prod(stepf%g2fft(g)))
               enddo
            else
               DO igptm = 1, mpdata%n_g(iq)
                  g = mpdata%g(:, mpdata%gptm_ptr(igptm, iq)) - g_t
                  cprod%data_c(hybdat%nbasp + igptm, iob + (iband - 1)*psize) = prod(stepf%g2fft(g))
               enddo
            endif
         enddo
      enddo
      !$OMP END DO
      deallocate (prod, psi_k)
      call fft%free()
      !$OMP END PARALLEL

      call stepf%free()

      call timestop("Big OMP loop")
      call psi_kqpt%free()
      call timestop("wavef_IS_FFT")
   end subroutine wavefproducts_IS_FFT

   subroutine wavef2rs(fi, lapw, zmat, gcutoff,  bandoi, bandof, jspin, psi)
!$    use omp_lib
      use m_types
      use m_fft_interface
      implicit none
      type(t_fleurinput), intent(in) :: fi
      type(t_lapw), intent(in)       :: lapw
      type(t_mat), intent(in)        :: zmat
      integer, intent(in)            :: jspin, bandoi, bandof
      real, intent(in)               :: gcutoff
      complex, intent(inout)         :: psi(0:, bandoi:) ! (nv,ne)

      type(t_fft) :: fft
      type(t_fftgrid) :: grid

      integer :: ivmap(SIZE(lapw%gvec, 2))
      integer :: iv, nu, n_threads

      psi = 0.0
      n_threads = 1
      !$OMP PARALLEL private(nu, iv, n_threads, fft, grid)  default(private) &
      !$OMP shared(fi, bandoi, bandof, zMat, psi,  ivmap, lapw, jspin, gcutoff)

      call grid%init(fi%cell, fi%sym, gcutoff)
      call fft%init(grid%dimensions, .false.)

      !$OMP DO
      do nu = bandoi, bandof
         call grid%putStateOnGrid(lapw, jspin, zMat, nu)
         psi(:,nu) = grid%grid
         call fft%exec(psi(:, nu))
      enddo
      !$OMP ENDDO
      call grid%free()
      call fft%free()
      !$OMP END PARALLEL
   end subroutine wavef2rs

   subroutine prep_list_of_gvec(lapw, mpdata, g_bounds, g_t, iq, jsp, pointer, gpt0, ngpt0)
      use m_types
      use m_juDFT
      implicit none
      type(t_lapw), intent(in)    :: lapw
      TYPE(t_mpdata), intent(in)         :: mpdata
      integer, intent(in)    :: g_bounds(:), g_t(:), iq, jsp
      integer, allocatable, intent(inout) :: pointer(:, :, :), gpt0(:, :)
      integer, intent(inout) :: ngpt0

      integer :: ic, ig1, igptm, iigptm, ok, g(3)

      allocate (pointer(-g_bounds(1):g_bounds(1), &
                        -g_bounds(2):g_bounds(2), &
                        -g_bounds(3):g_bounds(3)), stat=ok)
      IF (ok /= 0) call juDFT_error('wavefproducts_noinv2: error allocation pointer')
      allocate (gpt0(3, size(pointer)), stat=ok)
      IF (ok /= 0) call juDFT_error('wavefproducts_noinv2: error allocation gpt0')

      call timestart("prep list of Gvec")
      pointer = 0
      ic = 0
      DO ig1 = 1, lapw%nv(jsp)
         DO igptm = 1, mpdata%n_g(iq)
            iigptm = mpdata%gptm_ptr(igptm, iq)
            g = lapw%gvec(:, ig1, jsp) + mpdata%g(:, iigptm) - g_t
            IF (pointer(g(1), g(2), g(3)) == 0) THEN
               ic = ic + 1
               gpt0(:, ic) = g
               pointer(g(1), g(2), g(3)) = ic
            END IF
         END DO
      END DO
      ngpt0 = ic
      call timestop("prep list of Gvec")
   end subroutine prep_list_of_gvec

   function calc_number_of_basis_functions(lapw, atoms, noco) result(nbasfcn)
      use m_types
      implicit NONE
      type(t_lapw), intent(in)  :: lapw
      type(t_atoms), intent(in) :: atoms
      type(t_noco), intent(in)  :: noco
      integer                   :: nbasfcn

      if (noco%l_noco) then
         nbasfcn = lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot
      else
         nbasfcn = lapw%nv(1) + atoms%nlotot
      endif
   end function calc_number_of_basis_functions

   function outer_prod(x, y) result(outer)
      implicit NONE
      complex, intent(in) :: x(:), y(:)
      complex :: outer(size(x), size(y))
      integer  :: i, j

      do j = 1, size(y)
         do i = 1, size(x)
            outer(i, j) = x(i)*y(j)
         enddo
      enddo
   end function outer_prod
end module m_wavefproducts_aux
