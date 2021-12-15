module m_wavefproducts_inv
   USE m_types_hybdat
   use m_wavefproducts_noinv
   USE m_constants
   USE m_judft
   USE m_types
   USE m_types_hybinp
   USE m_util
   USE m_io_hybrid
   USE m_wrapper
   USE m_constants
   USE m_wavefproducts_aux

CONTAINS
   SUBROUTINE wavefproducts_inv(fi, ik, z_k, iq, jsp, bandoi, bandof, lapw, hybdat, mpdata, nococonv, stars, ikqpt, cmt_nk, cprod)
      IMPLICIT NONE
      type(t_fleurinput), intent(in):: fi
      TYPE(t_mpdata), intent(in)    :: mpdata
      type(t_nococonv), intent(in)  :: nococonv
      type(t_stars), intent(in)     :: stars
      type(t_mat), intent(in)       :: z_k  ! = z_k_p since ik < nkpt
      TYPE(t_lapw), INTENT(IN)      :: lapw
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      type(t_mat), intent(inout)    :: cprod

      ! - scalars -
      INTEGER, INTENT(IN)      :: jsp, ik, iq, bandoi, bandof
      INTEGER, INTENT(INOUT)   :: ikqpt
      complex, intent(in)  :: cmt_nk(:,:,:)


      ! - local scalars -
      INTEGER                 ::    g_t(3), i, j
      REAL                    ::    kqpt(3), kqpthlp(3)

      type(t_mat) ::  z_kqpt_p
      complex, allocatable :: c_phase_kqpt(:), tmp(:,:)

      CALL timestart("wavefproducts_inv")
      ikqpt = -1
      kqpthlp = fi%kpts%bkf(:, ik) + fi%kpts%bkf(:, iq)
      ! kqpt can lie outside the first BZ, transfer it back
      kqpt = fi%kpts%to_first_bz(kqpthlp)
      g_t = nint(kqpt - kqpthlp)

      ! determine number of kqpt
      ikqpt = fi%kpts%get_nk(kqpt)
      allocate (c_phase_kqpt(hybdat%nbands(fi%kpts%bkp(ikqpt),jsp)))
      IF (.not. fi%kpts%is_kpt(kqpt)) call juDFT_error('wavefproducts_inv5: k-point not found')

      !$acc data copyin(hybdat, hybdat%nbasp)
         !$acc data copyin(cprod) create(cprod%data_c) copyout(cprod%data_r)
            !$acc kernels present(cprod, cprod%data_r)
            cprod%data_r(:,:) = 0.0
            !$acc end kernels
            call wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                       ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)
         !$acc end data ! cprod

         
         allocate(tmp(hybdat%nbasp, cprod%matsize2))
         !$acc data copyout(tmp)
            !$acc kernels present(tmp)
            tmp = cmplx_0
            !$acc end kernels
            call wavefproducts_noinv_MT(fi, ik, iq, bandoi, bandof, nococonv, mpdata, hybdat, &
                                       jsp, ikqpt, z_kqpt_p, c_phase_kqpt, cmt_nk, tmp)
            call transform_to_realsph(fi, mpdata, tmp)
         !$acc end data
            
         call timestart("cpu cmplx2real copy")
         !$omp parallel do default(none) collapse(2) private(i,j) shared(cprod, hybdat, tmp)
         do i = 1,cprod%matsize2
            do j = 1,hybdat%nbasp
               cprod%data_r(j,i) = real(tmp(j,i))
            enddo 
         enddo
         !$omp end parallel do
         call timestop("cpu cmplx2real copy")
         deallocate(tmp)
      !$acc end data ! hybdat
      CALL timestop("wavefproducts_inv")
   END SUBROUTINE wavefproducts_inv

   subroutine transform_to_realsph(fi, mpdata, cprod)
      use m_constants
      implicit none 
      type(t_fleurinput), intent(in):: fi
      TYPE(t_mpdata), intent(in)    :: mpdata
      complex, intent(inout)        :: cprod(:,:)

      integer :: lm_0, lm, iatm, itype, l, m, partner

      !$acc data copyin(mpdata, mpdata%num_radbasfn)
         lm_0 = 0
         do iatm = 1,fi%atoms%nat 
            itype = fi%atoms%itype(iatm)

            ! The default(shared) in the OMP part of the following loop is needed to avoid compilation issues on gfortran 7.5.
            DO l = 0, fi%hybinp%lcutm1(itype)
               DO m = -l, l
                  lm = lm_0 + (m + l)*mpdata%num_radbasfn(l, itype)
                  if(m == 0) then
                     if(mod(l,2) == 1) then
                        !$acc kernels present(cprod, mpdata, mpdata%num_radbasfn)
                        cprod(lm+1:lm+mpdata%num_radbasfn(l, itype), :) &
                           = -ImagUnit * cprod(lm+1:lm+mpdata%num_radbasfn(l, itype), :)
                        !$acc end kernels
                     endif
                  else
                     if(m < 0) then
                        partner = lm + 2*abs(m)*mpdata%num_radbasfn(l,itype)
                        !$acc kernels present(cprod, mpdata, mpdata%num_radbasfn)
                        cprod(lm+1:lm+mpdata%num_radbasfn(l, itype), :) &
                           = (-1.0)**l * sqrt_2 * (-1.0)**m  * real(cprod(partner+1:partner+mpdata%num_radbasfn(l, itype), :))
                        !$acc end kernels
                     else 
                        !$acc kernels present(cprod, mpdata, mpdata%num_radbasfn)
                        cprod(lm+1:lm+mpdata%num_radbasfn(l, itype), :) &
                           = -(-1.0)**l * sqrt_2 * (-1.0)**m  * aimag(cprod(lm+1:lm+mpdata%num_radbasfn(l, itype), :))
                        !$acc end kernels
                     endif 
                  endif
               enddo 
               lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
            enddo
         enddo
      !$acc end data
   end subroutine transform_to_realsph
end module m_wavefproducts_inv