module m_wavefproducts_aux
   use m_types_fftGrid
   use m_types
CONTAINS
   subroutine wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                   ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)
      !$ use omp_lib
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
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat
      type(t_stars), intent(in)       :: stars
      type(t_mat), intent(in)         :: z_k
      type(t_mat), intent(inout)      :: z_kqpt_p, cprod
      !     - scalars -
      INTEGER, INTENT(IN)      ::  ik, iq, jsp, g_t(3), bandoi, bandof
      INTEGER, INTENT(IN)      ::  ikqpt
      !     - arrays -
      complex, intent(inout)    :: c_phase_kqpt(hybdat%nbands(ikqpt,jsp))

      complex, allocatable  :: prod(:,:), psi_k(:, :), psi_kqpt(:,:)

      type(t_mat)     :: z_kqpt
      type(t_lapw)    :: lapw_ikqpt
      type(t_fft)     :: fft, wavef2rs_fft
      type(t_fftgrid) :: stepf, grid


      integer, parameter :: blocksize = 512
      integer :: g(3), igptm, iob, n_omp, j, jstart, loop_length
      integer :: ok, nbasfcn, psize, iband, ierr, i, max_igptm
      integer, allocatable :: band_list(:), g_ptr(:)
      real    :: inv_vol, gcutoff, max_imag

      logical :: real_warned

      real_warned = .False.

      call timestart("wavef_IS_FFT")
      max_igptm = mpdata%n_g(iq)

      gcutoff = (2*fi%input%rkmax + fi%mpinp%g_cutoff) * fi%hybinp%fftcut
      inv_vol = 1/sqrt(fi%cell%omtil)
      psize = bandof - bandoi + 1
      !this is for the exact result. Christoph recommend 2*gmax+gcutm for later
      if (2*fi%input%rkmax + fi%mpinp%g_cutoff > fi%input%gmax) then
         write (*, *) "WARNING: not accurate enough: 2*kmax+gcutm >= fi%input%gmax"
         !call juDFT_error("not accurate enough: 2*kmax+gcutm >= fi%input%gmax")
      endif

      call stepf%init(fi%cell, fi%sym, gcutoff)
      block
         type(t_cell)         :: cell !unused 
         call stepf%putfieldOnGrid(stars, stars%ustep)
      end block
      call fft%init(stepf%dimensions, .false., batch_size=1, l_gpu=.True.)
      !$acc data copyin(stepf, stepf%grid, stepf%gridlength)
         ! after we transform psi_k*stepf*psi_kqpt back  to 
         ! G-space we have to divide by stepf%gridLength. We do this now

         !$acc kernels default(none) present(stepf, stepf%grid, stepf%gridLength)
         stepf%grid = stepf%grid * inv_vol / stepf%gridLength
         !$acc end kernels

         call fft%exec(stepf%grid)
         call fft%free()
         
         call setup_g_ptr(mpdata, stepf, g_t, iq, g_ptr)
         
         CALL lapw_ikqpt%init(fi, nococonv, ikqpt)

         nbasfcn = lapw_ikqpt%hyb_num_bas_fun(fi)
         call z_kqpt%alloc(z_k%l_real, nbasfcn, psize)
         call z_kqpt_p%init(z_kqpt)

         band_list = [(i, i=bandoi, bandof)]
         call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, fi%input, ikqpt, jsp, z_kqpt, &
                     c_phase=c_phase_kqpt, parent_z=z_kqpt_p, list=band_list)
#ifdef CPP_MPI
         call timestart("read_z barrier")
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         hybdat%max_q = hybdat%max_q - 1
         call timestop("read_z barrier")
#endif

         allocate(psi_kqpt(0:stepf%gridLength-1, psize), stat=ierr)
         if(ierr /= 0) call juDFT_error("can't alloc psi_kqpt")

         !$acc data create(psi_kqpt)
            call grid%init(fi%cell, fi%sym, gcutoff)
            call wavef2rs_fft%init(grid%dimensions, .false., batch_size=psize, l_gpu=.True.)
            !$acc data copyin(z_kqpt, z_kqpt%l_real, z_kqpt%data_r, z_kqpt%data_c, lapw_ikqpt, lapw_ikqpt%nv, lapw_ikqpt%gvec,&
            !$acc             jsp, bandoi, bandof, psize, grid, grid%dimensions)
               call timestart("1st wavef2rs")
               call wavef2rs(fi, lapw_ikqpt, z_kqpt, gcutoff, 1, psize, jsp, grid, wavef2rs_fft, psi_kqpt)
               call timestop("1st wavef2rs")

               !$acc kernels default(none) present(psi_kqpt, stepf, stepf%grid)
               do iob = 1, psize 
                  psi_kqpt(:,iob) = psi_kqpt(:,iob) * stepf%grid
               enddo
               !$acc end kernels
            !$acc end data
            call wavef2rs_fft%free()
            !call grid%free()

            call timestart("Big OMP loop")
#ifndef _OPENACC
!            !$OMP PARALLEL default(none) & ! Issue #770
!            !$OMP private(iband, iob, g, igptm, prod, psi_k, ok, fft, wavef2rs_fft, max_imag, grid) &
!            !$OMP shared(hybdat, psi_kqpt, cprod,  mpdata, iq, g_t, psize, gcutoff, max_igptm)&
!            !$OMP shared(jsp, z_k, stars, lapw, fi, inv_vol, ik, real_warned, n_omp, bandoi, stepf, g_ptr)
#endif

!            call timestart("alloc&init")
            allocate (prod(0:stepf%gridLength - 1, psize), stat=ok)
            if (ok /= 0) call juDFT_error("can't alloc prod")
            allocate (psi_k(0:stepf%gridLength - 1, 1), stat=ok)
            if (ok /= 0) call juDFT_error("can't alloc psi_k")

            call fft%init(stepf%dimensions, .true., batch_size=psize, l_gpu=.True.)
            call grid%init(fi%cell, fi%sym, gcutoff)
            call wavef2rs_fft%init(grid%dimensions, .false., batch_size=1, l_gpu=.True.)
!            call timestop("alloc&init")

            !$acc data copyin(z_k, z_k%l_real, z_k%data_r, z_k%data_c, lapw, lapw%nv, lapw%gvec)&
            !$acc      copyin(hybdat, hybdat%nbasp, g_ptr, grid, grid%dimensions, jsp)&
            !$acc      create(psi_k, prod)
#ifndef _OPENACC
!               !$OMP DO
#endif
               do iband = 1, hybdat%nbands(ik,jsp)
                  call wavef2rs(fi, lapw, z_k, gcutoff, iband, iband, jsp, grid, wavef2rs_fft, psi_k)
                  
                  !$acc kernels default(none) present(prod, psi_k, psi_kqpt, stepf, stepf%gridlength)               
                  do iob = 1, psize
                     do j = 0, stepf%gridlength-1
                        prod(j,iob) = conjg(psi_k(j, 1)) * psi_kqpt(j, iob)
                     enddo
                  enddo
                  !$acc end kernels

                  call fft%exec_batch(prod)
            
                  if (cprod%l_real) then
                     if (.not. real_warned) then
                        !$acc kernels present(prod) copyout(max_imag)
                        max_imag = maxval(abs(aimag(prod)))
                        !$acc end kernels
                        if(max_imag > 1e-8) then
                           write (*, *) "Imag part non-zero in too large"
                           real_warned = .True.
                        endif
                     endif
                        
                     !$acc kernels default(none) present(cprod, cprod%data_r, prod, g_ptr)
                     !$acc loop independent
                     do iob = 1, psize
                        !$acc loop independent
                        DO igptm = 1, max_igptm
                           cprod%data_r(hybdat%nbasp + igptm, iob + (iband - 1)*psize) = real(prod(g_ptr(igptm), iob))
                        enddo
                     enddo
                     !$acc end kernels
                  else
                     !$acc kernels default(none) present(cprod, cprod%data_c, prod, g_ptr)
                     !$acc loop independent
                     do iob = 1, psize
                        !$acc loop independent
                        DO igptm = 1, max_igptm
                           cprod%data_c(hybdat%nbasp + igptm, iob + (iband - 1)*psize) = prod(g_ptr(igptm), iob)
                        enddo
                     enddo
                     !$acc end kernels
                  endif
               enddo
#ifndef _OPENACC
!               !$OMP END DO
#endif
            !$acc end data 
            call fft%free()
            !call grid%free()
            call wavef2rs_fft%free()
         !$acc end data ! psi_kqpt
         deallocate (prod, psi_k)
      !$acc end data ! stepf, stepf%grid

#ifndef _OPENACC
!      !$OMP END PARALLEL
#endif
      !call stepf%free() 

      call timestop("Big OMP loop")
      deallocate(psi_kqpt)
      call timestop("wavef_IS_FFT")
   end subroutine wavefproducts_IS_FFT

   subroutine setup_g_ptr(mpdata, stepf, g_t, iq, g_out)
      implicit none
      type(t_mpdata), intent(in)          :: mpdata 
      type(t_fftgrid), intent(in)         :: stepf 
      integer, intent(in)                 :: g_t(:), iq
      integer, allocatable, intent(inout) :: g_out(:)

      integer :: igptm, g(3)

      if(allocated(g_out)) deallocate(g_out)
      allocate(g_out(mpdata%n_g(iq)))

      DO igptm = 1, mpdata%n_g(iq)
         g = mpdata%g(:, mpdata%gptm_ptr(igptm, iq)) - g_t
         g_out(igptm) = stepf%g2fft(g)
      enddo
   end subroutine setup_g_ptr

   subroutine wavef2rs(fi, lapw, zmat, gcutoff,  bandoi, bandof, jspin, grid, fft, psi)
      ! put block of wave functions through FFT
!$    use omp_lib
      use m_types
      use m_fft_interface
      implicit none
      type(t_fleurinput), intent(in) :: fi
      type(t_lapw), intent(in)       :: lapw
      type(t_mat), intent(in)        :: zmat
      integer, intent(in)            :: jspin, bandoi, bandof
      real, intent(in)               :: gcutoff
      type(t_fftgrid), intent(inout) :: grid
      type(t_fft), intent(inout)     :: fft
      complex, intent(inout)         :: psi(0:, bandoi:) ! (nv,ne)

      integer :: iv, nu, psize, dims(3)

#ifndef _OPENACC 
      !$omp parallel do default(none) private(nu) shared(grid, bandoi, bandof, lapw, jspin, zMat, psi)
#endif
      do nu = bandoi, bandof
         call grid%put_state_on_external_grid(lapw, jspin, zMat, nu, psi(:,nu), l_gpu=.True.)
      enddo
#ifndef _OPENACC 
      !$omp end parallel do
#endif   

      call fft%exec_batch(psi)
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
