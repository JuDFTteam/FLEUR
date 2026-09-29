module m_wavefproducts_inv
   USE m_vac_rows, ONLY: NVAC_MPB
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
      type(t_mat)  :: z_kq
      type(t_lapw) :: lapw_kq
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

      !$acc data copyin(hybdat, hybdat%n_mt)
         !$acc data copyin(cprod) create(cprod%data_c) copyout(cprod%data_r)
            !$acc kernels present(cprod, cprod%data_r)
            cprod%data_r(:,:) = 0.0
            !$acc end kernels
            if (fi%input%film) then
               call wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                          ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod, &
                                          lapw_kq_out=lapw_kq, z_kq_out=z_kq)
            else
               call wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                          ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)
            endif
         !$acc end data ! cprod

         if (fi%input%film) then
            block
               use m_wavefproducts_vac, only: wavefproducts_vac
               complex, allocatable :: tmp_vac(:,:)
               ! films: vacuum rows, computed complex and rotated to the real basis
               allocate(tmp_vac(cprod%matsize1, cprod%matsize2), source=cmplx_0)
               call wavefproducts_vac(fi, ik, iq, ikqpt, g_t, jsp, bandoi, bandof, mpdata, &
                                      hybdat, lapw, lapw_kq, z_k, z_kq, tmp_vac)
               call vac_to_realbasis(fi, mpdata, hybdat, iq, tmp_vac, cprod)
               deallocate(tmp_vac)
            end block
         endif

         
         allocate(tmp(hybdat%n_mt, cprod%matsize2))
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
            do j = 1,hybdat%n_mt
               cprod%data_r(j,i) = real(tmp(j,i))
            enddo 
         enddo
         !$omp end parallel do
         call timestop("cpu cmplx2real copy")
         deallocate(tmp)
      !$acc end data ! hybdat

      CALL timestop("wavefproducts_inv")
   END SUBROUTINE wavefproducts_inv

   !>Vacuum rows of cprod in the real pair (M_1 +- M_2)/sqrt(2); with inversion symmetry
   !>vac_abcof gives conjugate coefficients in the two vacua, so both combinations are real.
   subroutine vac_to_realbasis(fi, mpdata, hybdat, iq, cprod_c, cprod)
      use m_vac_rows, only: row_offset, basfn_offset
      use m_constants, only: sqrt_2
      implicit none
      type(t_fleurinput), intent(in) :: fi
      type(t_mpdata), intent(in)     :: mpdata
      type(t_hybdat), intent(in)     :: hybdat
      integer, intent(in)            :: iq
      complex, intent(in)            :: cprod_c(:,:)
      type(t_mat), intent(inout)     :: cprod

      integer :: igm, ilen, nn, b1, b2, i, j

      if (NVAC_MPB /= 2) call juDFT_error( &
         "invs film: the real vacuum basis needs both vacua (nvac = 2)", &
         calledby="wavefproducts_inv", hint="nvac = 1 is not covered yet")

      do igm = 1, mpdata%n_g_vac(iq)
         ilen = mpdata%glen_ptr_vac(igm, iq)
         nn = mpdata%num_zbasfn_vac(ilen, 1)
         if (nn == 0) cycle
         b1 = row_offset(mpdata, hybdat, iq, 1) + basfn_offset(mpdata, iq, 1, igm)
         b2 = row_offset(mpdata, hybdat, iq, 2) + basfn_offset(mpdata, iq, 2, igm)
         !$omp parallel do default(none) collapse(2) private(i,j) &
         !$omp shared(cprod, cprod_c, b1, b2, nn)
         do i = 1, cprod%matsize2
            do j = 1, nn
               cprod%data_r(b1 + j, i) = sqrt_2*real(cprod_c(b1 + j, i))
               cprod%data_r(b2 + j, i) = sqrt_2*aimag(cprod_c(b1 + j, i))
            enddo
         enddo
         !$omp end parallel do
      enddo
   end subroutine vac_to_realbasis
   subroutine transform_to_realsph(fi, mpdata, cprod)
      use m_constants
      implicit none 
      type(t_fleurinput), intent(in):: fi
      TYPE(t_mpdata), intent(in)    :: mpdata
      complex, intent(inout)        :: cprod(:,:)

      integer :: lm_0, lm, iatm, iatm2, itype, l, m, partner, ioffset, ishift

      !$acc data copyin(mpdata, mpdata%num_radbasfn)
         lm_0 = 0
         do iatm = 1,fi%atoms%nat 
            itype = fi%atoms%itype(iatm)

            iatm2 = fi%sym%invsatnr(iatm)
            IF(iatm2.EQ.0) iatm2 = iatm

            ioffset = sum((/((2*l + 1)*mpdata%num_radbasfn(l, itype), l=0, fi%hybinp%lcutm1(itype))/))
            
            IF(iatm2.LT.iatm) THEN ! iatm is the second of two atoms that are mapped onto each other by inversion symmetry
               DO l = 0, fi%hybinp%lcutm1(itype)
                  lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
               END DO
               CYCLE
            ELSE IF (iatm2.GT.iatm) THEN
               ! In this case we already make everything correct in wavefproducts_noinv_MT
               DO l = 0, fi%hybinp%lcutm1(itype)
                  lm = lm_0
                  DO m = -l, l

                     ishift = -2 * m * mpdata%num_radbasfn(l, itype)
!                     lm1 = lm + (iatm - fi%atoms%firstAtom(itype))*ioffset
!                     lm2 = lm + (iatm2 - fi%atoms%firstAtom(itype))*ioffset + ishift

!                     DO i = 1, mpdata%num_radbasfn(l, itype)
!                        cprod(i + lm1, (j-bandoi+1) + (k-1)*psize) = REAL(cprod(i + lm1, (j-bandoi+1) + (k-1)*psize))
!                     END DO

                     lm = lm + mpdata%num_radbasfn(l, itype)
                  END DO
                  lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
               END DO
            ELSE
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
            END IF
         enddo
      !$acc end data
   end subroutine transform_to_realsph
end module m_wavefproducts_inv
