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
   use m_matmul_dgemm
   LOGICAL, PARAMETER:: zero_order = .false., ibs_corr = .false.

CONTAINS
   SUBROUTINE exchange_valence_hf(k_pack, fi, fmpi, z_k, mpdata, jsp, hybdat, lapw, eig_irr, &
                                  n_q, wl_iks, xcpot, nococonv, stars, nsest, indx_sest, cmt_nk, mat_ex)
      
      USE m_wrapper
      USE m_trafo
      USE m_wavefproducts
      USE m_olap
      USE m_hsefunctional
      USE m_io_hybrid
      USE m_kp_perturbation
      use m_spmm_inv
      use m_spmm_noinv
      use m_work_package
      use m_judft 
#ifdef CPP_MPI
      use mpi
#endif
#ifdef _OPENACC
      USE cublas
#define CPP_zgemm cublaszgemm
#define CPP_dgemm cublasdgemm
#else
#define CPP_zgemm zgemm
#define CPP_dgemm dgemm
#endif
      IMPLICIT NONE

      type(t_k_package), intent(in)     :: k_pack
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      type(t_mat), intent(in)           :: z_k
      TYPE(t_xcpot_inbuild), INTENT(IN) :: xcpot
      TYPE(t_mpdata), intent(inout)     :: mpdata
      TYPE(t_nococonv), INTENT(IN)      :: nococonv
      TYPE(t_lapw), INTENT(IN)          :: lapw
      type(t_stars), intent(in)         :: stars
      TYPE(t_mat), INTENT(INOUT)        :: mat_ex
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      ! scalars
      INTEGER, INTENT(IN)    :: jsp

      ! arrays
      INTEGER, INTENT(IN)    ::  n_q(:)
      INTEGER, INTENT(IN)    ::  nsest(:)
      INTEGER, INTENT(IN)    ::  indx_sest(:, :)

      complex, intent(in)    :: cmt_nk(:,:,:)

      REAL, INTENT(IN)    ::  eig_irr(:, :)
      REAL, INTENT(IN)    ::  wl_iks(:, :)

      ! local scalars
      INTEGER                 ::  iband, iband1, jq, iq, nq_idx
      INTEGER                 ::  i, ierr, ik
      INTEGER                 ::  j, iq_p, start, stride
      INTEGER                 ::  n1, n2, nn2, me, max_band_pack
      INTEGER                 ::  ikqpt, iob, m, n, k, lda, ldb, ldc
      INTEGER                 ::  ok, psize, n_parts, ipart, ibando

      REAL, SAVE             ::  divergence

      COMPLEX                 ::  cdum2
      COMPLEX                 ::  exch0

      LOGICAL, SAVE           ::  initialize = .true.

      ! local arrays
      COMPLEX                          :: exchcorrect(fi%kpts%nkptf)
      COMPLEX, allocatable             :: dcprod(:,:,:) ! (hybdat%nbands(k_pack%nk,jsp), hybdat%nbands(k_pack%nk,jsp), 3)
      COMPLEX, allocatable             :: exch_vv(:,:) !(hybdat%nbands(k_pack%nk,jsp), hybdat%nbands(k_pack%nk,jsp))
      COMPLEX                          :: hessian(3, 3)
      COMPLEX                          :: proj_ibsc(3, MAXVAL(hybdat%nobd(:, jsp)), hybdat%nbands(k_pack%nk,jsp))
      COMPLEX                          :: olap_ibsc(3, 3, MAXVAL(hybdat%nobd(:, jsp)), MAXVAL(hybdat%nobd(:, jsp)))
      COMPLEX, ALLOCATABLE  :: phase_vv(:, :), c_coul_wavf(:,:), dot_result_c(:,:)
      REAL, ALLOCATABLE     :: r_coul_wavf(:,:), dot_result_r(:,:)
      LOGICAL                          :: occup(fi%input%neig), conjg_mtir

#define CPP_cprod_r cprod_vv%data_r
#define CPP_cprod_c cprod_vv%data_c

      type(t_mat)          :: cprod_vv, carr3_vv
      CALL timestart("valence exchange calculation")
      ik = k_pack%nk

      IF (initialize) THEN !it .eq. 1 .and. ik .eq. 1) THEN
         call calc_divergence(fi%cell, fi%kpts, divergence)
         if(fmpi%irank == 0) write (*,*) "Divergence:", divergence
         initialize = .false.
      END IF

      ! calculate valence-valence-valence-valence, core-valence-valence-valence
      ! and core-valence-valence-core exchange at current k-point
      ! the sum over the inner occupied valence states is restricted to the EIBZ(k)
      ! the contribution of the Gamma-point is treated separately (see below)

      call timestart("alloc phase_vv & dot_res")
      if(mat_ex%l_real) then
         allocate(dot_result_r(hybdat%nbands(ik,jsp), hybdat%nbands(ik,jsp)), stat=ierr, source=0.0)
         allocate (phase_vv(MAXVAL(hybdat%nobd(:, jsp)), hybdat%nbands(ik,jsp)), stat=ok, source=cmplx_0)
      else 
         allocate(dot_result_c(hybdat%nbands(ik,jsp), hybdat%nbands(ik,jsp)), stat=ierr, source=cmplx_0)
         allocate (phase_vv(0,0), stat=ok, source=cmplx_0)
      endif
      IF (ok /= 0) call judft_error('exchange_val_hf: error allocation phase')
      if(ierr /= 0) call judft_error("can't alloc dot_res")

      allocate(exch_vv(hybdat%nbands(ik,jsp), hybdat%nbands(ik,jsp)), stat=ierr, source=cmplx_0)
      if(ierr /= 0) call judft_error("Can't alloc exch_vv")


      !$acc data copyout(exch_vv) copyin(hybdat, hybdat%nbands, hybdat%nbasm, nsest, indx_sest) 
         !$acc kernels present(exch_vv) default(none)
         exch_vv = 0
         !$acc end kernels

         call timestop("alloc phase_vv & dot_res")
         
         call timestart("q_loop")

         DO jq = 1, size(k_pack%q_packs)
            call timestart("initial setup")
            iq = k_pack%q_packs(jq)%ptr
            iq_p = fi%kpts%bkp(iq)

            ikqpt = fi%kpts%get_nk(fi%kpts%to_first_bz(fi%kpts%bkf(:, ik) + fi%kpts%bkf(:, iq)))

            n_parts = size(k_pack%q_packs(jq)%band_packs)
            start  = k_pack%q_packs(jq)%submpi%rank + 1
            stride = k_pack%q_packs(jq)%submpi%size
            call timestop("initial setup")
            do ipart = start, n_parts, stride
               !if (n_parts > 1) write (*, *) "Part ("//int2str(ipart)//"/"//int2str(n_parts)//") ik= "//int2str(ik)//" jq= "//int2str(jq)
               psize = k_pack%q_packs(jq)%band_packs(ipart)%psize
               ibando = k_pack%q_packs(jq)%band_packs(ipart)%start_idx
               call timestart("alloc & copy")
               call cprod_vv%alloc(mat_ex%l_real, hybdat%nbasm(iq), psize*hybdat%nbands(ik,jsp), mat_name="cprod_vv")
               call alloc_dev_cpy(cprod_vv, CPP_cprod_r, CPP_cprod_c)
               call timestop("alloc & copy")
               IF (mat_ex%l_real) THEN
                  CALL wavefproducts_inv(fi, ik, z_k, iq, jsp, ibando, ibando + psize - 1, lapw, &
                                       hybdat, mpdata, nococonv, stars, ikqpt, cmt_nk, cprod_vv)
               ELSE
                  CALL wavefproducts_noinv(fi, ik, z_k, iq, jsp, ibando, ibando + psize - 1, lapw,&
                                          hybdat, mpdata, nococonv, stars, ikqpt, cmt_nk, cprod_vv)
               END IF
               ! The sparse matrix technique is not feasible for the HSE
               ! functional. Thus, a dynamic adjustment is implemented
               ! The mixed basis functions and the potential difference
               ! are Fourier transformed, so that the exchange can be calculated
               ! in Fourier space
               !! REIMPLEMENTING (notes in lab book)
               IF (xcpot%is_name("hse") .OR. xcpot%is_name("vhse")) THEN
                  CALL timestart("hse: dynamic hse adjustment")
                  iband1 = hybdat%nobd(ikqpt, jsp)
                  exch_vv = exch_vv + &
                            dynamic_hse_adjustment(fi%atoms, fi%kpts%bkf(:, iq), iq, &
                                                   fi%kpts%nkptf, fi%cell%bmat, fi%cell%omtil, &
                                                   fi%hybinp%lcutm1, maxval(fi%hybinp%lcutm1), mpdata%num_radbasfn, maxval(mpdata%num_radbasfn), mpdata%g, &
                                                   mpdata%n_g(iq), mpdata%gptm_ptr(:, iq), mpdata%num_gpts(), mpdata%radbasfn_mt, &
                                                   hybdat%nbasm(iq), iband1, hybdat%nbands(ik,jsp), nsest, ibando, psize, indx_sest, &
                                                   fi%sym, fmpi%irank, cprod_vv%data_r(:, :), &
                                                   cprod_vv%data_c(:, :), mat_ex%l_real, wl_iks(:iband1, ikqpt), n_q(jq))
                  CALL timestop("hse: dynamic hse adjustment")
               END IF

               ! the Coulomb matrix is only evaluated at the irrecuible k-points
               ! bra_trafo transforms cprod instead of rotating the Coulomb matrix
               ! from IBZ to current k-point

               call timestart("bra_trafo stuff")
               IF (fi%kpts%bkp(iq) /= iq) THEN
                  call carr3_vv%init(cprod_vv, mat_name="carr_3")
                  call bra_trafo(fi, mpdata, hybdat, hybdat%nbands(ik,jsp), iq, psize, phase_vv, cprod_vv, carr3_vv)
                  call cprod_vv%copy(carr3_vv, 1, 1)
                  call carr3_vv%free()
               ELSE
                  phase_vv(:, :) = cmplx_1
               END IF
               call timestop("bra_trafo stuff")
               
               call timestart("alloc coul_wavf")
               if(cprod_vv%l_real) then 
                  allocate(r_coul_wavf(cprod_vv%matsize1, cprod_vv%matsize2), stat=ierr)
                  allocate(c_coul_wavf(0,0))
               else
                  allocate(c_coul_wavf(cprod_vv%matsize1, cprod_vv%matsize2), stat=ierr)
                  allocate(r_coul_wavf(0,0))
               endif
               if(ierr /= 0) call judft_error("can't alloc coul_wavf")
               r_coul_wavf = 0.0
               call timestop("alloc coul_wavf")
               
               call timestart("exchange matrix")
               call timestart("sparse matrix products")
               IF (mat_ex%l_real) THEN
                  call spmm_invs(fi, mpdata, hybdat, iq_p, cprod_vv%data_r, r_coul_wavf)
               ELSE
                  conjg_mtir = (fi%kpts%bksym(iq) > fi%sym%nop)
                  call spmm_noinvs(fi, mpdata, hybdat, iq_p, conjg_mtir, cprod_vv%data_c, c_coul_wavf)
               END IF
               call timestop("sparse matrix products")

               !$acc enter data copyin(phase_vv, r_coul_wavf, c_coul_wavf)
               nq_idx = k_pack%q_packs(jq)%rank


               call timestart("apply prefactors carr1_v")
               !$acc data copyin(psize, wl_iks, n_q, nq_idx, ibando, ikqpt)
                  if (mat_ex%l_real) then
#ifdef _OPENACC
                     call timestart("cpy cprod")
                     call dlacpy("N", size(cprod_vv%data_r, 1), size(cprod_vv%data_r, 2), cprod_vv%data_r, size(cprod_vv%data_r, 1), CPP_cprod_r, size(CPP_cprod_r,1))
                     call timestop("cpy cprod")
#endif
                     !$acc enter data copyin(CPP_cprod_r)

                     !$acc parallel loop default(none) collapse(3) private(iband, iob, i)&
                     !$acc present(r_coul_wavf, hybdat, hybdat%nbands, hybdat%nbasm, psize, wl_iks, ikqpt, ibando)&
                     !$acc present(n_q, nq_idx)
                     DO iband = 1, hybdat%nbands(ik,jsp)
                        DO iob = 1, psize
                           do i = 1, hybdat%nbasm(iq)
                              r_coul_wavf(i, iob + psize*(iband - 1)) = r_coul_wavf(i, iob + psize*(iband - 1)) * wl_iks(ibando + iob - 1, ikqpt) / n_q(nq_idx)                        
                           enddo
                        enddo
                     enddo
                     !$acc end parallel loop

                  else
#ifdef _OPENACC
                     call timestart("cpy cprod")
                     call zlacpy("N", size(cprod_vv%data_c, 1), size(cprod_vv%data_c, 2), cprod_vv%data_c, size(cprod_vv%data_c, 1), CPP_cprod_c, size(CPP_cprod_c,1))
                     call timestop("cpy cprod")
#endif
                     !$acc enter data copyin(CPP_cprod_c)

                     !$acc parallel loop default(none) collapse(3) private(iband, iob, i)&
                     !$acc present(c_coul_wavf, hybdat, hybdat%nbands, hybdat%nbasm, psize, wl_iks, ikqpt, ibando)&
                     !$acc present(n_q, nq_idx)
                     DO iband = 1, hybdat%nbands(ik,jsp)
                        DO iob = 1, psize
                           do i = 1, hybdat%nbasm(iq)
                              c_coul_wavf(i, iob + psize*(iband - 1))  = c_coul_wavf(i, iob + psize*(iband - 1)) * wl_iks(ibando + iob - 1, ikqpt)/n_q(nq_idx)
                           enddo
                        enddo
                     enddo
                     !$acc end parallel loop
                  endif
               !$acc end data ! (psize, wl_iks, n_q, nq_idx, ibando, ikqpt)
               call timestop("apply prefactors carr1_v")

               call timestart("exch_vv dot prod")
               m = hybdat%nbands(ik,jsp)
               n = hybdat%nbands(ik,jsp)
               k = hybdat%nbasm(iq)
               lda = hybdat%nbasm(iq)*psize
               ldb = hybdat%nbasm(iq)*psize
               ldc = hybdat%nbands(ik,jsp)

               IF (mat_ex%l_real) THEN
                  !calculate all dotproducts for the current iob -> need to skip intermediate iob
                  !$acc enter data create(dot_result_r) 
                  DO iob = 1, psize
                     call timestart("CPP_dgemm")
                     !call blas_matmul(m,n,k,r_coul_wavf(:,iob:),CPP_cprod_r(:, iob:),dot_result_r,op_a="T")
                     ASSOCIATE(prod_data=>CPP_cprod_r) !persuade NVHPC that it knows CPP_CPROD_r
                        !$acc host_data use_device(r_coul_wavf, prod_data, dot_result_r)
                        call CPP_dgemm("T", "N", m, n, k, 1.0, r_coul_wavf(1, iob), lda, prod_data(1, iob), ldb, 0.0, dot_result_r , ldc)
                        !$acc end host_data
                     end ASSOCIATE   
                     !$acc wait
                     call timestop("CPP_dgemm")

                     !$acc kernels present(exch_vv, dot_result_r, phase_vv, hybdat, hybdat%nbands, nsest, indx_sest) default(none)
                     DO iband = 1, hybdat%nbands(ik,jsp)
                        DO n2 = 1, nsest(iband)
                           nn2 = indx_sest(n2, iband)
                           exch_vv(nn2, iband) = exch_vv(nn2, iband) + phase_vv(iob, nn2)*conjg(phase_vv(iob, iband))*dot_result_r(iband, nn2)
                        enddo
                     END DO
                     !$acc end kernels
                  END DO
                  !$acc exit data delete(CPP_cprod_r, dot_result_r)
               ELSE
                  !calculate all dotproducts for the current iob -> need to skip intermediate iob
                  !$acc enter data create(dot_result_c) 
                  DO iob = 1, psize
                     call timestart("CPP_zgemm")
                     ASSOCIATE(prod_data=>CPP_cprod_c) !persuade NVHPC that it knows CPP_CPROD_C
                        !$acc host_data use_device(c_coul_wavf, prod_data, dot_result_c)
                        call CPP_zgemm("C", "N", m, n, k, cmplx_1, c_coul_wavf(1, iob), lda, prod_data(1, iob), ldb, cmplx_0, dot_result_c, ldc)
                        !$acc end host_data
                     end ASSOCIATE   
                     !$acc wait
                     call timestop("CPP_zgemm")

                     !$acc kernels present(exch_vv, dot_result_c, hybdat, hybdat%nbands, nsest, indx_sest)  default(none)
                     DO iband = 1, hybdat%nbands(ik,jsp)
                        DO n2 = 1, nsest(iband)
                           nn2 = indx_sest(n2, iband)
                           exch_vv(nn2, iband) = exch_vv(nn2, iband) + dot_result_c(iband, nn2)
                        enddo
                     END DO
                     !$acc end kernels
                  enddo
                  !$acc exit data delete(CPP_cprod_c, dot_result_c)
               END IF
               !$acc exit data delete(phase_vv, r_coul_wavf, c_coul_wavf)
               call timestop("exch_vv dot prod")
               call timestop("exchange matrix")
               
               call cprod_vv%free()
               if(allocated(r_coul_wavf)) deallocate(r_coul_wavf)
               if(allocated(c_coul_wavf)) deallocate(c_coul_wavf)
            enddo
         END DO  !jq
      !$acc end data ! exch_vv hybdat, hybdat%nbands, hybdat%nbasm, nsest, indx_sest
      call timestop("q_loop")

      if(allocated(dot_result_r)) deallocate(dot_result_r)
      if(allocated(dot_result_c)) deallocate(dot_result_c)

      ! add contribution of the gamma point to the different cases (exch_vv,exch_cv,exch_cc)

      ! valence-valence-valence-valence exchange
      call timestart("gamma point treatment")
      IF ((.not. xcpot%is_name("hse")) .AND. &
          (.not. xcpot%is_name("vhse")) .AND. &
          k_pack%submpi%root()) THEN ! no gamma point correction needed for HSE functional

         IF (zero_order .and. .not. ibs_corr) THEN
            WRITE (oUnit, '(A)') ' Take zero order terms into account.'
         ELSE IF (zero_order .and. ibs_corr) THEN
            WRITE (oUnit, '(A)') ' Take zero order terms and ibs-correction into account.'
         END IF

         IF (zero_order) THEN
            CALL dwavefproducts(dcprod, ik, 1, hybdat%nbands(ik,jsp), 1, hybdat%nbands(ik,jsp), .false., fi%input, fi%atoms, mpdata, fi%hybinp, &
                                fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, lapw,  jsp, eig_irr)

            ! make dcprod hermitian
            DO n1 = 1, hybdat%nbands(ik,jsp)
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
         DO i = 1, hybdat%results%neig(ik, jsp)
            IF (hybdat%results%ef >= eig_irr(i, ik)) THEN
               occup(i) = .true.
            ELSE IF ((eig_irr(i, ik) - hybdat%results%ef) <= 1E-06) THEN
               occup(i) = .true.
            END IF
         END DO

         DO n1 = 1, hybdat%nbands(ik,jsp)
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
                        DO iband = 1, hybdat%nbands(ik,jsp)
                           IF (occup(iband)) THEN
                              hessian(i, j) = hessian(i, j) + conjg(dcprod(iband, n1, i))*dcprod(iband, nn2, j)
                           END IF
                           hessian(i, j) = hessian(i, j) - dcprod(iband, nn2, i)*conjg(dcprod(iband, n1, j))
                        END DO

                        ! ibs correction
                        IF (ibs_corr) THEN
                           hessian(i, j) = hessian(i, j) - olap_ibsc(i, j, n1, nn2)/fi%cell%omtil
                           DO iband = 1, hybdat%nbands(ik,jsp)
                              hessian(i, j) = hessian(i, j) + conjg(proj_ibsc(i, nn2, iband))*proj_ibsc(j, n1, iband)/fi%cell%omtil
                           END DO
                        END IF
                     END DO
                  ELSE
                     DO i = 1, 3
                        j = i
                        DO iband = 1, hybdat%nbands(ik,jsp)
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
      call timestop("gamma point treatment")

      IF (mat_ex%l_real) THEN
         IF (any(abs(aimag(exch_vv)) > 1E-08)) then 
            CALL judft_warn('unusually large imaginary part of exch_vv. Max:' // float2str(maxval(abs(aimag(exch_vv)))), &
                                                               calledby='exchange_val_hf.F90')
         endif
      END IF

      ! write exch_vv in mat_ex
      call timestart("alloc mat_ex")
      if (.not. mat_ex%allocated()) then
         if (k_pack%submpi%root()) then
            CALL mat_ex%alloc(matsize1=hybdat%nbands(ik,jsp))
         else
            CALL mat_ex%alloc(matsize1=1)
         endif
      endif
      call timestop("alloc mat_ex")

#ifdef CPP_MPI 
      call timestart("pre exchmat reduce barrier")
      call MPI_Barrier(k_pack%submpi%comm, ierr)
      call timestop("pre exchmat reduce barrier")
#endif

      call timestart("reduce exch_vv>mat_ex")
      IF (mat_ex%l_real) THEN
#ifdef CPP_MPI
         call MPI_Reduce(real(exch_vv), mat_ex%data_r, hybdat%nbands(ik,jsp)**2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, k_pack%submpi%comm, ierr)
#else
         mat_ex%data_r = exch_vv
#endif
      ELSE
#ifdef CPP_MPI
         call MPI_Reduce(exch_vv, mat_ex%data_c, hybdat%nbands(ik,jsp)**2, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, k_pack%submpi%comm, ierr)
#else
         mat_ex%data_c = exch_vv
#endif
      END IF
      call timestop("reduce exch_vv>mat_ex")
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

      call timestart("calc_divergence")
      expo = 5e-3
      rrad = sqrt(-log(5e-3)/expo)
      cdum = sqrt(expo)*rrad
      divergence = real(cell%omtil/(tpi_const**2)*sqrt(pi_const/expo)*cerf(cdum))
      rrad = rrad**2
      nkpt3 = kpts%calcNkpt3()
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
      call timestop("calc_divergence")
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
      REAL    :: rrad, k(3), knorm2
      COMPLEX :: cdum

      rrad = sqrt(-log(5e-3)/expo)
      cdum = sqrt(expo)*rrad
      divergence = real(cell%omtil/(tpi_const**2)*sqrt(pi_const/expo)*cerf(cdum))
      rrad = rrad**2

      do ikpt = 1, kpts%nkptf
         k = kpts%bkf(:, ikpt)
         knorm2 = norm2(k)
         IF (knorm2 < rrad) THEN
            divergence = divergence - exp(-expo*knorm2)/knorm2/kpts%nkptf
         END IF
      enddo
   end function calc_divergence2

   subroutine recombine_parts(in_part, ipart, psizes, out_total)
      use m_types
      type(t_mat), intent(in)    :: in_part
      integer, intent(in)        :: ipart, psizes(:)
      type(t_mat), intent(inout) :: out_total

      integer :: nbands, iband, iob, offset, i, tsize
      logical :: l_real

      l_real = in_part%l_real
      tsize = sum(psizes)

      nbands = in_part%matsize2/psizes(ipart)
      if (out_total%matsize2/tsize /= nbands) call judft_error("nbands seems different")
      offset = 0
      do i = 1, ipart - 1
         offset = offset + psizes(i)
      enddo

      do iband = 1, nbands
         do iob = 1, psizes(ipart)
            if (l_real) then
               out_total%data_r(:, iob + (iband - 1)*tsize + offset) = in_part%data_r(:, iob + (iband - 1)*psizes(ipart))
            else
               out_total%data_c(:, iob + (iband - 1)*tsize + offset) = in_part%data_c(:, iob + (iband - 1)*psizes(ipart))
            endif
         enddo
      enddo
   end subroutine recombine_parts

   subroutine alloc_dev_cpy(mat, arr_r, arr_c)
      implicit none 
      type(t_mat), intent(in)             :: mat

      real, allocatable, intent(inout)    :: arr_r(:,:)
      complex, allocatable, intent(inout) :: arr_c(:,:) 
      integer :: ierr

#ifdef _OPENACC
      if(allocated(arr_r)) deallocate(arr_r)
      if(allocated(arr_c)) deallocate(arr_c)

      if(mat%l_real) then 
         allocate(arr_r(mat%matsize1, mat%matsize2), stat=ierr)
      else
         allocate(arr_c(mat%matsize1, mat%matsize2), stat=ierr)
      endif 

      if(ierr /= 0) call judft_error("can't alloc. prob. no mem.")
#endif
   end subroutine alloc_dev_cpy

END MODULE m_exchange_valence_hf
