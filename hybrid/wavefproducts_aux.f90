module m_wavefproducts_aux

CONTAINS
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

   subroutine wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
      ikqpt, z_k, c_phase_k, z_kqpt_p, c_phase_kqpt, cprod)
      use m_types
      use m_constants
      use m_judft
      use m_fft_interface
      use m_io_hybinp
      use m_juDFT
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
      complex, intent(inout)    :: c_phase_k(hybdat%nbands(ik)), c_phase_kqpt(hybdat%nbands(ikqpt))
      
      complex  :: prod(0:27*stars%mx1*stars%mx2*stars%mx3-1)
      
      type(t_mat)               :: z_kqpt
      type(t_lapw)              :: lapw_ikqpt
      integer :: length_zfft(3), g(3), igptm, gshift(3), iob
      integer :: ok, ne, nbasfcn, fftd, psize, iband, irs, ob
      integer, allocatable :: iob_arr(:), iband_arr(:)
      real    :: q(3), inv_vol 
      type(t_mat)  :: psi_k, psi_kqpt

      call timestart("wavef_IS_FFT")
      inv_vol = 1/sqrt(fi%cell%omtil)
      length_zfft = [3*stars%mx1, 3*stars%mx2, 3*stars%mx3]
      fftd = product(length_zfft)
      psize = bandof - bandoi + 1
      !this is for the exact result. Christoph recommend 2*gmax+gcutm for later
      if(2*fi%input%rkmax + fi%mpinp%g_cutoff > fi%input%gmax) then 
         write (*,*) "WARNING: not accurate enough: 2*kmax+gcutm >= fi%input%gmax"
         !call juDFT_error("not accurate enough: 2*kmax+gcutm >= fi%input%gmax")
      endif
      
      CALL lapw_ikqpt%init(fi, nococonv, ikqpt)
      nbasfcn = lapw_ikqpt%hyb_num_bas_fun(fi)
      call z_kqpt%alloc(z_k%l_real, nbasfcn, fi%input%neig)
      call z_kqpt_p%init(z_kqpt)
      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, fi%input, ikqpt, jsp, z_kqpt, &
                  c_phase=c_phase_kqpt, parent_z=z_kqpt_p)
      
      call psi_k%alloc(.false., fftd, hybdat%nbands(ik))
      call psi_kqpt%alloc(.false., fftd, psize)

      call wavef2rs(fi, lapw, stars, z_k, length_zfft, 1, hybdat%nbands(ik), jsp, psi_k%data_c)
      call wavef2rs(fi, lapw_ikqpt, stars, z_kqpt, length_zfft, bandoi, bandof, jsp, psi_kqpt%data_c)
      
      !$OMP PARALLEL DO default(none) &
      !$OMP private(iband, irs) shared(psi_k, stars, inv_vol)&
      !$OMP collapse(2)
      do iband = 1, hybdat%nbands(ik)
         do irs = 1,fftd
            psi_k%data_c(irs,iband) = conjg(psi_k%data_c(irs,iband)) * stars%ufft(irs-1) * inv_vol
         enddo 
      enddo
      !$OMP END PARALLEL DO


      !$OMP PARALLEL DO default(none) &
      !$OMP private(iband, iob, g, igptm, prod) &
      !$OMP shared(hybdat, psi_k, psi_kqpt, cprod, length_zfft, mpdata, iq, g_t, psize) &
      !$OMP collapse(2)
      do iband = 1,hybdat%nbands(ik)
         do iob = 1, psize 
            prod = psi_k%data_c(:,iband) * psi_kqpt%data_c(:,iob)
            call fft_interface(3, length_zfft, prod, .true.)
            if(cprod%l_real) then
               if(any(abs(aimag(prod)) > 1e-10)) call juDFT_error("Imag part non-zero in is_fft")
            endif
            
            ! we still have to devide by the number of mesh points
            prod = prod / product(length_zfft)

            if(cprod%l_real) then
               DO igptm = 1, mpdata%n_g(iq)
                  g = mpdata%g(:, mpdata%gptm_ptr(igptm, iq)) - g_t
                  cprod%data_r(hybdat%nbasp+igptm, iob + (iband-1)*psize) = real(prod(g2fft(length_zfft,g)))        
               enddo
            else
               DO igptm = 1, mpdata%n_g(iq)
                  g = mpdata%g(:, mpdata%gptm_ptr(igptm, iq)) - g_t
                  cprod%data_c(hybdat%nbasp+igptm, iob + (iband-1)*psize) = prod(g2fft(length_zfft,g))        
               enddo
            endif  
         enddo 
      enddo
      !$OMP END PARALLEL DO

      !$OMP barrier
      call timestop("wavef_IS_FFT")
   end subroutine wavefproducts_IS_FFT

   subroutine wavef2rs(fi, lapw, stars, zmat, length_zfft, bandoi, bandof, jspin, psi)
      use m_types
      use m_fft_interface
      implicit none
      type(t_fleurinput), intent(in) :: fi
      type(t_lapw), intent(in)       :: lapw
      type(t_mat), intent(in)        :: zmat
      type(t_stars), intent(in)      :: stars
      integer, intent(in)            :: jspin, bandoi, bandof, length_zfft(3)
      complex, intent(inout)         :: psi(0:,bandoi:) ! (nv,ne)

      integer :: ivmap(SIZE(lapw%gvec, 2))
      integer :: iv, nu

      DO iv = 1, lapw%nv(jspin)
         ivmap(iv) = g2fft(length_zfft, lapw%gvec(:, iv, jspin))
      ENDDO

      psi = 0.0
      do nu = bandoi, bandof
         !------> map WF into FFTbox
         DO iv = 1, lapw%nv(jspin)
            if(zMat%l_real) then
               psi(ivmap(iv), nu) = zMat%data_r(iv, nu)
            else
               psi(ivmap(iv), nu) = zMat%data_c(iv, nu)
            endif
         ENDDO

         call fft_interface(3, length_zfft, psi(:,nu), .false., ivmap(1:lapw%nv(jspin)))
      enddo
   end subroutine wavef2rs

end module m_wavefproducts_aux
