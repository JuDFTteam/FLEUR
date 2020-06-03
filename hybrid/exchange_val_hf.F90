!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!     Calculates the HF exchange term
!
!                                          s          s*          s            s*
!                                       phi    (r) phi     (r) phi     (r') phi    (r')
!                         occ.             n_1k       n'k+q       n'k+q        n_2k
!     exchange(n,q)  =  - SUM  INT INT  ------------------------------------------- dr dr'
!                         k,n'                           | r - r' |
!
!                         occ                  s          s    ~        ~       s         s
!                    =  - SUM  SUM  v     < phi      | phi     M    > < M    phi     | phi      >
!                         k,n' I,J   k,IJ      n'k+q      n_1k  q,I      q,J    n_2k      n'k+q
!
!     for the different combinations of n_1 and n_2 and where n' runs only over the valence states.
!     ( n_1,n_2:  valence-valence, core-core,core-valence )
!
!
!     At the Gamma point (k=0) v diverges. After diagonalization of v at k=0 the divergence is
!     restricted to the head element I=1. Furthermore, we expand <...> with kp perturbation theory.
!     As a result, the total I=1 element is given by a sum of a divergent 1/k**2-term and an
!     angular dependent term. The former is separated from the numerical k-summation and treated
!     analytically while the latter is spherically averaged and added to the k=0 contribution of
!     the numerical k-summation. (A better knowledge of the integrand's behavior at the BZ edges
!     might further improve the integration.)
!
!     The divergence at the Gamma point is integrated with one of the following algorithms:
! (1) Switching-Off Function
!     In a sphere of radius k0=radshmin/2 a switching-off function g(k)=1-(k/k0)**n*(n+1-n*k/k0)
!     (n=npot) is defined. The 1/k**2 divergence is subtracted from the BZ integral in the form
!     g(k)/k**2 and integrated analytically. The non-divergent rest is integrated numerically.
! (2) Periodic Function (similar to the one used by Massidda PRB 48, 5058)
!     The function  F(k) = SUM(G) exp(-expo*|k+G|**3) / |k+G|**2  is subtracted from the BZ integral
!     and integrated analytically. The non-divergent rest is integrated numerically.
!     The parameter expo is chosen such that exp(-expo*q**3)=1/2
!     with q = radius of sphere with same volume as BZ.
! (3) Periodic Function (same as Massidda's) with expo->0
!     The function  F(k) = lim(expo->0) SUM(G) exp(-expo*|k+G|**2) / |k+G|**2  is subtracted from
!     the BZ integral and integrated analytically. The contribution to the BZ integral including
!     the "tail" is
!     vol/(8*pi**3) INT F(k) d^3k - P SUM(k) F(k)  ( P = principal value ) .
!     For expo->0 the two terms diverge. Therefore a cutoff radius q0 is introduced and related to
!     expo by exp(-expo*q0**2)=delta  ( delta = small value, e.g., delta = 1*10.0**-10 ) .
!     The resulting formula
!     vol/(4*pi**1.5*sqrt(expo)) * erf(sqrt(a)*q0) - sum(q,0<q<q0) exp(-expo*q**2)/q**2
!     converges well with q0. (Should be the default.)

MODULE m_exchange_valence_hf

   USE m_constants
   USE m_types
   USE m_util

   LOGICAL, PARAMETER:: zero_order = .false., ibs_corr = .false.
   INTEGER, PARAMETER:: maxmem = 600

CONTAINS
   SUBROUTINE exchange_valence_hf(ik, fi, z_k, c_phase_k, mpdata, jsp, hybdat, lapw, eig_irr, results, &
                                  pointer_EIBZ, n_q, wl_iks, xcpot, nococonv, stars, nsest, indx_sest, mpi, mat_ex)

      USE m_wrapper
      USE m_trafo
      USE m_wavefproducts
      USE m_olap
      USE m_spmvec
      USE m_hsefunctional
      USE m_io_hybinp
      USE m_kp_perturbation
      use m_spmm
      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      type(t_mat), intent(in)           :: z_k
      TYPE(t_results), INTENT(IN)       :: results
      TYPE(t_xcpot_inbuild), INTENT(IN) :: xcpot
      TYPE(t_mpi), INTENT(IN)           :: mpi
      TYPE(t_mpdata), intent(inout)     :: mpdata
      TYPE(t_nococonv), INTENT(IN)      :: nococonv
      TYPE(t_lapw), INTENT(IN)          :: lapw
      type(t_stars), intent(in)         :: stars
      TYPE(t_mat), INTENT(INOUT)        :: mat_ex
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      ! blas 
      real, external      :: ddot
      complex, external   :: zdotc

      ! scalars
      INTEGER, INTENT(IN)    :: jsp
      INTEGER, INTENT(IN)    :: ik

      ! arrays
      INTEGER, INTENT(IN)    ::  n_q(:)

      INTEGER, INTENT(IN)    ::  pointer_EIBZ(:)
      INTEGER, INTENT(IN)    ::  nsest(:)
      INTEGER, INTENT(IN)    ::  indx_sest(:, :)

      REAL, INTENT(IN)    ::  eig_irr(:, :)
      REAL, INTENT(IN)    ::  wl_iks(:, :)
      complex, intent(in) :: c_phase_k(hybdat%nbands(ik))

      ! local scalars
      INTEGER                 ::  iband, iband1, jq, iq
      INTEGER                 ::  i, ierr
      INTEGER                 ::  j, iq_p
      INTEGER                 ::  n1, n2, nn2
      INTEGER                 ::  ikqpt, iob, m,n,k,lda,ldb,ldc
      INTEGER                 ::  ok, psize, n_parts, ipart, ibando
      integer, allocatable    :: start_idx(:), psizes(:)

      REAL, SAVE             ::  divergence

      COMPLEX                 ::  cdum, cdum2
      COMPLEX                 ::  exch0

      LOGICAL, SAVE           ::  initialize = .true.

      ! local arrays
      COMPLEX              :: exchcorrect(fi%kpts%nkptf)
      COMPLEX              :: dcprod(hybdat%nbands(ik), hybdat%nbands(ik), 3)
      COMPLEX              :: exch_vv(hybdat%nbands(ik), hybdat%nbands(ik))
      COMPLEX              :: hessian(3, 3), ctmp
      COMPLEX              :: proj_ibsc(3, MAXVAL(hybdat%nobd(:, jsp)), hybdat%nbands(ik))
      COMPLEX              :: olap_ibsc(3, 3, MAXVAL(hybdat%nobd(:, jsp)), MAXVAL(hybdat%nobd(:, jsp)))
      COMPLEX, ALLOCATABLE :: phase_vv(:, :)
      REAL                 :: kqpt(3), kqpthlp(3), target_psize, rtmp

      LOGICAL              :: occup(fi%input%neig), conjg_mtir
      type(t_mat)          :: carr1_v, cprod_vv, carr3_vv, dot_result
      character(len=300)   :: errmsg
      CALL timestart("valence exchange calculation")

      IF (initialize) THEN !it .eq. 1 .and. ik .eq. 1) THEN
         call calc_divergence(fi%cell, fi%kpts, divergence)
         PRINT *, "Divergence:", divergence
         initialize = .false.
      END IF

      ! calculate valence-valence-valence-valence, core-valence-valence-valence
      ! and core-valence-valence-core exchange at current k-point
      ! the sum over the inner occupied valence states is restricted to the EIBZ(k)
      ! the contribution of the Gamma-point is treated separately (see below)

      allocate (phase_vv(MAXVAL(hybdat%nobd(:, jsp)), hybdat%nbands(ik)), stat=ok, source=cmplx_0)
      call dot_result%alloc(mat_ex%l_real, hybdat%nbands(ik), hybdat%nbands(ik))
      IF (ok /= 0) call judft_error('exchange_val_hf: error allocation phase')

      exch_vv = 0

      DO jq = 1,fi%kpts%EIBZ(ik)%nkpt
         iq = pointer_EIBZ(jq)
         iq_p = fi%kpts%bkp(iq)


         ikqpt = fi%kpts%get_nk(fi%kpts%to_first_bz(fi%kpts%bkf(:,ik) + fi%kpts%bkf(:,iq)))
         ! arrays should be less than 5 gb
         if(mat_ex%l_real) then
            target_psize = 5e9/( 8.0 * maxval(hybdat%nbasm) * hybdat%nbands(ik)) 
         else
            target_psize = 5e9/(16.0 * maxval(hybdat%nbasm) * hybdat%nbands(ik)) 
         endif
         n_parts = ceiling(hybdat%nobd(ikqpt, jsp)/target_psize)
         call split_iob_loop(hybdat, hybdat%nobd(ikqpt, jsp), n_parts, start_idx, psizes)
         do ipart = 1, n_parts
            if(n_parts > 1) write (*,*) "Part (" // int2str(ipart) //"/"// int2str(n_parts) // ")"
            psize = psizes(ipart)
            ibando = start_idx(ipart)
            call cprod_vv%alloc(mat_ex%l_real, hybdat%nbasm(iq), psize * hybdat%nbands(ik))

            IF (mat_ex%l_real) THEN
               CALL wavefproducts_inv(fi, ik, z_k, iq, jsp, ibando, ibando+psize-1, lapw, hybdat, mpdata, nococonv, stars, ikqpt, cprod_vv)
            ELSE
               CALL wavefproducts_noinv(fi, ik, z_k, iq, jsp, ibando, ibando+psize-1, lapw, hybdat, mpdata, nococonv, stars, ikqpt, cprod_vv)
            END IF

            ! The sparse matrix technique is not feasible for the HSE
            ! functional. Thus, a dynamic adjustment is implemented
            ! The mixed basis functions and the potential difference
            ! are Fourier transformed, so that the exchange can be calculated
            ! in Fourier space
            IF (xcpot%is_name("hse") .OR. xcpot%is_name("vhse")) THEN
               call judft_error("HSE not implemented")
               ! iband1 = hybdat%nobd(ikqpt, jsp)

               ! exch_vv = exch_vv + &
               !           dynamic_hse_adjustment(fi%atoms%rmsh, fi%atoms%rmt, fi%atoms%dx, fi%atoms%jri, fi%atoms%jmtd, fi%kpts%bkf(:, iq), iq, &
               !                                  fi%kpts%nkptf, fi%cell%bmat, fi%cell%omtil, fi%atoms%ntype, fi%atoms%neq, fi%atoms%nat, fi%atoms%taual, &
               !                                  fi%hybinp%lcutm1, maxval(fi%hybinp%lcutm1), mpdata%num_radbasfn, maxval(mpdata%num_radbasfn), mpdata%g, &
               !                                  mpdata%n_g(iq), mpdata%gptm_ptr(:, iq), mpdata%num_gpts(), mpdata%radbasfn_mt, &
               !                                  hybdat%nbasm(iq), iband1, hybdat%nbands(ik), nsest, 1, MAXVAL(hybdat%nobd(:, jsp)), indx_sest, &
               !                                  fi%sym%invsat, fi%sym%invsatnr, mpi%irank, cprod_vv_r(:hybdat%nbasm(iq), :, :), &
               !                                  cprod_vv_c(:hybdat%nbasm(iq), :, :), mat_ex%l_real, wl_iks(:iband1, ikqpt), n_q(jq))
            END IF

            ! the Coulomb matrix is only evaluated at the irrecuible k-points
            ! bra_trafo transforms cprod instead of rotating the Coulomb matrix
            ! from IBZ to current k-point
            IF (fi%kpts%bkp(iq) /= iq) THEN
               call carr3_vv%init(cprod_vv)
               call bra_trafo(fi, mpdata, hybdat, hybdat%nbands(ik), iq, jsp, psize, phase_vv, cprod_vv, carr3_vv)
               call cprod_vv%copy(carr3_vv, 1,1)
               call carr3_vv%free()
            ELSE
               phase_vv(:, :) = cmplx_1
            END IF

            call carr1_v%init(cprod_vv)
            ! calculate exchange matrix at iq
            call timestart("exchange matrix")
            ! finish coulomb bcast
            call hybdat%coul(iq_p)%mpi_wait()
            call timestart("sparse matrix products")
            IF (mat_ex%l_real) THEN
               call spmm_invs(fi, mpdata, hybdat, iq_p, cprod_vv, carr1_v)
            ELSE
               conjg_mtir = (fi%kpts%bksym(iq) > fi%sym%nop)
               call spmm_noinvs(fi, mpdata, hybdat, iq_p, conjg_mtir, cprod_vv, carr1_v)
            END IF
            call timestop("sparse matrix products")

            DO iband = 1, hybdat%nbands(ik)
               call timestart("apply prefactors carr1_v")
               if(mat_ex%l_real) then
                  DO iob = 1, psize
                     do i=1,hybdat%nbasm(iq)
                        carr1_v%data_r(i,iob + psize*(iband-1)) = carr1_v%data_r(i,iob + psize*(iband-1)) * wl_iks(ibando+iob-1, ikqpt) * conjg(phase_vv(iob, iband))/n_q(jq)
                     enddo
                  enddo
               else
                  DO iob = 1, psize
                     do i=1,hybdat%nbasm(iq)
                        carr1_v%data_c(i,iob + psize*(iband-1)) = carr1_v%data_c(i,iob + psize*(iband-1)) * wl_iks(ibando+iob-1, ikqpt) * conjg(phase_vv(iob, iband))/n_q(jq)
                     enddo
                  enddo
               endif
               call timestop("apply prefactors carr1_v")
            enddo

            call timestart("exch_vv dot prod")
            m = hybdat%nbands(ik)
            n = hybdat%nbands(ik)
            k = hybdat%nbasm(iq)
            lda = hybdat%nbasm(iq) * psize
            ldb = hybdat%nbasm(iq) * psize
            ldc = hybdat%nbands(ik)
            IF (mat_ex%l_real) THEN
               !calculate all dotproducts for the current iob -> need to skip intermediate iob
               DO iob = 1, psize
                  call dgemm("T", "N", m, n, k, 1.0, carr1_v%data_r(1, iob), lda, cprod_vv%data_r(1, iob), ldb, 0.0, dot_result%data_r, ldc)

                  DO iband = 1, hybdat%nbands(ik)
                     DO n2 = 1, nsest(iband)
                        nn2 = indx_sest(n2, iband)
                        exch_vv(nn2, iband) = exch_vv(nn2, iband) + phase_vv(iob, nn2) * dot_result%data_r(iband, nn2)
                     enddo
                  END DO 
               END DO  
            ELSE
               !calculate all dotproducts for the current iob -> need to skip intermediate iob
               DO iob = 1, psize
                  call zgemm("C", "N", m, n, k, cmplx_1, carr1_v%data_c(1, iob), lda, cprod_vv%data_c(1, iob), ldb, cmplx_0, dot_result%data_c, ldc)

                  DO iband = 1, hybdat%nbands(ik)
                     DO n2 = 1, nsest(iband)
                        nn2 = indx_sest(n2, iband)
                        exch_vv(nn2, iband) = exch_vv(nn2, iband) + phase_vv(iob, nn2) * dot_result%data_c(iband, nn2)
                     enddo
                  END DO 
               enddo
            END IF
            call timestop("exch_vv dot prod")

            call timestop("exchange matrix")

            call cprod_vv%free()
            call carr1_v%free()
         enddo
      END DO  !jq
      call dot_result%free()

!   WRITE(7001,'(a,i7)') 'ik: ', ik
!   DO n1=1,hybdat%nbands(ik)
!      DO n2=1,n1
!         WRITE(7001,'(2i7,2f15.8)') n2, n1, exch_vv(n2,n1)
!     END DO
!   END DO

      ! add contribution of the gamma point to the different cases (exch_vv,exch_cv,exch_cc)

      ! valence-valence-valence-valence exchange

      IF ((.not. xcpot%is_name("hse")) .AND. (.not. xcpot%is_name("vhse"))) THEN ! no gamma point correction needed for HSE functional
         IF (zero_order .and. .not. ibs_corr) THEN
            WRITE (oUnit, '(A)') ' Take zero order terms into account.'
         ELSE IF (zero_order .and. ibs_corr) THEN
            WRITE (oUnit, '(A)') ' Take zero order terms and ibs-correction into account.'
         END IF

         IF (zero_order) THEN
            CALL dwavefproducts(dcprod, ik, 1, hybdat%nbands(ik), 1, hybdat%nbands(ik), .false., fi%input, fi%atoms, mpdata, fi%hybinp, &
                                fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, lapw, fi%oneD, jsp, eig_irr)

            ! make dcprod hermitian
            DO n1 = 1, hybdat%nbands(ik)
               DO n2 = 1, n1
                  dcprod(n1, n2, :) = (dcprod(n1, n2, :) - conjg(dcprod(n2, n1, :)))/2
                  dcprod(n2, n1, :) = -conjg(dcprod(n1, n2, :))
               END DO
            END DO

            IF (ibs_corr) THEN
               CALL ibs_correction(ik, fi%atoms, fi%input, jsp, hybdat, mpdata, fi%hybinp, &
                                   lapw, fi%kpts, fi%cell, MAXVAL(hybdat%nobd(:, jsp)), &
                                   fi%sym, fi%noco, nococonv, proj_ibsc, olap_ibsc)
            END IF
         END IF

         !This should be done with w_iks I guess!TODO
         occup = .false.
         DO i = 1, hybdat%ne_eig(ik)
            IF (results%ef >= eig_irr(i, ik)) THEN
               occup(i) = .true.
            ELSE IF ((eig_irr(i, ik) - results%ef) <= 1E-06) THEN
               occup(i) = .true.
            END IF
         END DO

         DO n1 = 1, hybdat%nbands(ik)
            DO n2 = 1, nsest(n1)!n1
               nn2 = indx_sest(n2, n1)
               exchcorrect = 0
               exch0 = 0

               ! if zero_order = .true. add averaged k-dependent term to the numerical integration at Gamma-point contribution

               ! if we start with a system with a small DFT band gap (like GaAs), the contribution
               ! of the highest occupied and lowest unoccupied state in Hessian is typically
               ! large; a correct numerical integration requires a dense k-point mesh, so
               ! we don't add the contribution exchcorrect for such materials

               IF (zero_order) THEN
                  hessian = 0
                  IF (occup(n1) .and. occup(nn2)) THEN
                     DO i = 1, 3
                        j = i
                        DO iband = 1, hybdat%nbands(ik)
                           IF (occup(iband)) THEN
                              hessian(i, j) = hessian(i, j) + conjg(dcprod(iband, n1, i))*dcprod(iband, nn2, j)
                           END IF
                           hessian(i, j) = hessian(i, j) - dcprod(iband, nn2, i)*conjg(dcprod(iband, n1, j))
                        END DO

                        ! ibs correction
                        IF (ibs_corr) THEN
                           hessian(i, j) = hessian(i, j) - olap_ibsc(i, j, n1, nn2)/fi%cell%omtil
                           DO iband = 1, hybdat%nbands(ik)
                              hessian(i, j) = hessian(i, j) + conjg(proj_ibsc(i, nn2, iband))*proj_ibsc(j, n1, iband)/fi%cell%omtil
                           END DO
                        END IF
                     END DO
                  ELSE
                     DO i = 1, 3
                        j = i
                        DO iband = 1, hybdat%nbands(ik)
                           IF (occup(iband)) THEN
                              hessian(i, j) = hessian(i, j) + conjg(dcprod(iband, n1, i))*dcprod(iband, nn2, j)
                           END IF
                        END DO
                     END DO
                  END IF

                  exchcorrect(1) = fpi_const/3*(hessian(1, 1) + hessian(2, 2) + hessian(3, 3))
                  exch0 = exchcorrect(1)/fi%kpts%nkptf
               END IF

               ! tail correction/contribution from all other k-points (it  goes into exchcorrect )

               ! Analytic contribution

               cdum2 = 0
               !multiply divergent contribution with occupation number;
               !this only affects metals
               IF (n1 == nn2) THEN
                  cdum2 = fpi_const/fi%cell%omtil*divergence*wl_iks(n1, ik)*fi%kpts%nkptf
               END IF

               ! due to the symmetrization afterwards the factor 1/n_q(1) must be added

               IF (n1 == nn2) hybdat%div_vv(n1, ik, jsp) = REAL(cdum2)
               exch_vv(nn2, n1) = exch_vv(nn2, n1) + (exch0 + cdum2)/n_q(1)

            END DO !n2
         END DO !n1
      END IF ! xcpot%icorr .ne. icorr_hse

      IF (mat_ex%l_real) THEN
         IF (any(abs(aimag(exch_vv)) > 1E-08)) CALL judft_warn('unusally large imaginary part of exch_vv', &
                                                               calledby='exchange_val_hf.F90')
      END IF

!   WRITE(7000,'(a,i7)') 'ik: ', ik
!   DO n1=1,hybdat%nbands(ik)
!      DO n2=1,n1
!         WRITE(7000,'(2i7,2f15.8)') n2, n1, exch_vv(n2,n1)
!      END DO
!   END DO

      ! write exch_vv in mat_ex
      CALL mat_ex%alloc(matsize1=hybdat%nbands(ik))
      IF (mat_ex%l_real) THEN
         mat_ex%data_r = exch_vv
      ELSE
         mat_ex%data_c = exch_vv
      END IF
      CALL timestop("valence exchange calculation")

   END SUBROUTINE exchange_valence_hf

   SUBROUTINE calc_divergence(cell, kpts, divergence)

      IMPLICIT NONE

      TYPE(t_cell), INTENT(IN)  :: cell
      TYPE(t_kpts), INTENT(IN)  :: kpts
      REAL, INTENT(INOUT) :: divergence

      INTEGER :: ix, iy, iz, sign, n
      logical :: found
      REAL    :: expo, rrad, k(3), kv1(3), kv2(3), kv3(3), knorm2, nkpt3(3)
      COMPLEX :: cdum

      expo = 5e-3
      rrad = sqrt(-log(5e-3)/expo)
      cdum = sqrt(expo)*rrad
      divergence = cell%omtil/(tpi_const**2)*sqrt(pi_const/expo)*cerf(cdum)
      rrad = rrad**2
      nkpt3 = kpts%nkpt3()
      kv1 = cell%bmat(1, :)/nkpt3(1)
      kv2 = cell%bmat(2, :)/nkpt3(2)
      kv3 = cell%bmat(3, :)/nkpt3(3)
      n = 1
      found = .true.

      DO WHILE (found)
         found = .false.
         DO ix = -n, n
            DO iy = -(n - abs(ix)), n - abs(ix)
               iz = n - abs(ix) - abs(iy)
               DO sign = -1, 1, 2
                  iz = sign*iz
                  k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
                  k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
                  k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
                  knorm2 = k(1)**2 + k(2)**2 + k(3)**2
                  IF (knorm2 < rrad) THEN
                     found = .true.
                     divergence = divergence - exp(-expo*knorm2)/knorm2/kpts%nkptf
                  END IF
                  IF (iz == 0) exit
               END DO
            END DO
         END DO
         n = n + 1
      END DO

   END SUBROUTINE calc_divergence

   function calc_divergence2(cell, kpts) result(divergence)
      USE m_types
      USE m_constants
      USE m_util, ONLY: cerf
      implicit none
      TYPE(t_cell), INTENT(IN)  :: cell
      TYPE(t_kpts), INTENT(IN)  :: kpts
      REAL                      :: divergence

      INTEGER :: ikpt
      REAL, PARAMETER :: expo = 5e-3
      REAL    :: rrad, k(3), kv1(3), kv2(3), kv3(3), knorm2
      COMPLEX :: cdum

      rrad = sqrt(-log(5e-3)/expo)
      cdum = sqrt(expo)*rrad
      divergence = cell%omtil/(tpi_const**2)*sqrt(pi_const/expo)*cerf(cdum)
      rrad = rrad**2

      do ikpt = 1, kpts%nkptf
         k = kpts%bkf(:, ikpt)
         knorm2 = norm2(k)
         IF (knorm2 < rrad) THEN
            divergence = divergence - exp(-expo*knorm2)/knorm2/kpts%nkptf
         END IF
      enddo
   end function calc_divergence2

   subroutine split_iob_loop(hybdat, n_total, n_parts, start_idx, psize)
      use m_types
      implicit none
      type(t_hybdat), intent(inout)       :: hybdat
      integer, intent(in)                 :: n_total, n_parts
      integer, allocatable, intent(inout) :: start_idx(:), psize(:)

      integer             :: n_loops, i, big_size, small_size, end_idx

      if(allocated(start_idx)) deallocate(start_idx)
      if(allocated(psize)) deallocate(psize)
      allocate(start_idx(n_parts), psize(n_parts))

      small_size = floor((1.0*n_total)/n_parts)
      big_size = small_size +1

      end_idx = 0
      do i = 1,n_parts
         psize(i) = merge(big_size, small_size,i <= mod(n_total, n_parts))

         start_idx(i) = end_idx + 1
         end_idx = start_idx(i) + psize(i) - 1
      enddo
      if(hybdat%l_print_iob_splitting) then
         write (*,*) "Split iob loop into " // int2str(n_parts) // " parts"
         write (*,*) "sizes: ", psize(1), psize(n_parts)
         hybdat%l_print_iob_splitting = .False.
      endif
   end subroutine split_iob_loop

   subroutine recombine_parts(in_part, ipart, psizes, out_total)
      use m_types 
      type(t_mat), intent(in)    :: in_part
      integer, intent(in)        :: ipart, psizes(:)
      type(t_mat), intent(inout) :: out_total 
      
      integer :: nbands, iband, iob, offset, i, tsize
      logical :: l_real 

      l_real = in_part%l_real
      tsize = sum(psizes)

      nbands = in_part%matsize2 / psizes(ipart)
      if(out_total%matsize2 / tsize /= nbands) call judft_error("nbands seems different")
      offset = 0 
      do i=1,ipart-1 
         offset = offset + psizes(i)
      enddo

      do iband = 1, nbands 
         do iob = 1,psizes(ipart)
            if(l_real) then 
               out_total%data_r(:,iob + (iband-1)*tsize + offset) = in_part%data_r(:,iob + (iband-1) * psizes(ipart))
            else
               out_total%data_c(:,iob + (iband-1)*tsize + offset) = in_part%data_c(:,iob + (iband-1) * psizes(ipart))
            endif 
         enddo
      enddo
   end subroutine recombine_parts
END MODULE m_exchange_valence_hf
