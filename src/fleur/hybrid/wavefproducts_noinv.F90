module m_wavefproducts_noinv
   USE m_types_hybdat

CONTAINS
   SUBROUTINE wavefproducts_noinv(fi, ik, z_k, iq, jsp, bandoi, bandof, lapw, hybdat, mpdata, nococonv, stars, ikqpt, cmt_nk, cprod)
      USE m_types
      use m_juDFT
      use m_wavefproducts_aux
      use m_constants, only: cmplx_0
      IMPLICIT NONE

      type(t_fleurinput), intent(in)  :: fi
      type(t_nococonv), intent(in)    :: nococonv
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_mpdata), intent(in)      :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat
      type(t_mat), intent(in)         :: z_k ! z_k is also z_k_p since ik < nkpt
      type(t_stars), intent(in)       :: stars
      type(t_mat), intent(inout)      :: cprod

!     - scalars -
      INTEGER, INTENT(IN)        ::  ik, iq, jsp, bandoi, bandof
      INTEGER, INTENT(INOUT)     ::  ikqpt

      complex, intent(in)  :: cmt_nk(:,:,:)

      INTEGER              :: g_t(3)
      REAL                 :: kqpt(3), kqpthlp(3)
      complex, allocatable :: c_phase_kqpt(:)
      type(t_mat)          :: z_kqpt_p

      call timestart("wavefproducts_noinv")
      ! calculate ikqpt
      kqpthlp = fi%kpts%bkf(:, ik) + fi%kpts%bkf(:, iq)
      kqpt = fi%kpts%to_first_bz(kqpthlp)

      ! if k+q outside of first BZ put we need this shift
      g_t = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      ikqpt = fi%kpts%get_nk(kqpt)
      call timestart("alloc c_phase_kqpt")
      allocate (c_phase_kqpt(hybdat%nbands(ikqpt,jsp)), source=cmplx_0)
      call timestop("alloc c_phase_kqpt")

      IF (.not. fi%kpts%is_kpt(kqpt)) then
         call juDFT_error('wavefproducts: k-point not found')
      endif

      !$acc data copyin(cprod) create(cprod%data_r) copyout(cprod%data_c)
         !$acc kernels 
         cprod%data_c(:,:) = 0.0
         !$acc end kernels
         call wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                    ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)

         call wavefproducts_noinv_MT(fi, ik, iq, bandoi, bandof, nococonv, mpdata, hybdat, &
                                    jsp, ikqpt, z_kqpt_p, c_phase_kqpt, cmt_nk, cprod%data_c)
      !$acc end data ! cprod

      call timestop("wavefproducts_noinv")

   END SUBROUTINE wavefproducts_noinv

   subroutine wavefproducts_noinv_MT(fi, ik, iq, bandoi, bandof, nococonv, mpdata, hybdat, jsp, ikqpt, &
                                     z_kqpt_p, c_phase_kqpt, cmt_nk, cprod)
      use m_types
      USE m_constants
      use m_io_hybrid
      use m_judft
      use m_wavefproducts_aux
      use m_calc_cmt
      IMPLICIT NONE
      type(t_fleurinput), intent(in)  :: fi
      type(t_nococonv), intent(in)    :: nococonv
      TYPE(t_mpdata), INTENT(IN)      :: mpdata
      TYPE(t_hybdat), INTENT(IN)      :: hybdat
      type(t_mat), intent(in)         :: z_kqpt_p
      complex, intent(inout)          :: cprod(:,:)

      !     - scalars -
      INTEGER, INTENT(IN)     ::  ik, iq, jsp, bandoi, bandof
      INTEGER, INTENT(IN)     ::  ikqpt

      !     - arrays -
      complex, intent(in)     :: c_phase_kqpt(hybdat%nbands(ikqpt,jsp))

      complex, intent(in)    :: cmt_nk(:,:,:)

      !     - local scalars -
      INTEGER                 ::  iatm, iatm2, l, n, l1, l2, n1, n2, lm_0, lm1_0, lm2_0, lm1_cprod, lm2_cprod
      INTEGER                 ::  lm, lm1, lm2, m1, m2, i, ll, j, k, ok
      INTEGER                 ::  itype, ieq, m, psize, ishift, ioffset

      REAL                    ::  sr2, atom_phase1, atom_phase2

      COMPLEX                 ::  atom_phase,  cscal, cscal2, add1, add2

      LOGICAL                 ::  offdiag

      !      - local arrays -
      INTEGER                 ::  lmstart(0:fi%atoms%lmaxd, fi%atoms%ntype)

      COMPLEX, allocatable    ::  cmt_ikqpt(:,:,:)

      call timestart("wavefproducts_noinv5 MT")
      allocate(cmt_ikqpt(bandoi:bandof, hybdat%maxlmindx, fi%atoms%nat), stat=ok, source=cmplx_0)
      if(ok /= 0) call juDFT_error("alloc cmt_ikqpt")
      
      psize = bandof-bandoi+1
      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      call timestart("set lmstart")
      DO itype = 1, fi%atoms%ntype
         DO l = 0, fi%atoms%lmax(itype)
            lmstart(l, itype) = sum([(mpdata%num_radfun_per_l(ll, itype)*(2*ll + 1), ll=0, l - 1)])
         END DO
      END DO
      call timestop("set lmstart")

      sr2 = SQRT(2.0)

      ! read in cmt coefficients from direct access file cmt
      call calc_cmt(fi%atoms, fi%cell, fi%input, fi%noco, nococonv, fi%hybinp, hybdat, mpdata, fi%kpts, fi%sym,  z_kqpt_p, jsp, ikqpt, c_phase_kqpt, cmt_ikqpt)

      call timestart("loop over l, l1, l2, n, n1, n2")
      !$acc data copyin(mpdata, mpdata%num_radbasfn, mpdata%num_radfun_per_l, mpdata%l1, mpdata%l2, mpdata%n1, mpdata%n2,&
      !$acc        hybdat, hybdat%prodm, hybdat%nbands, hybdat%nindxp1, hybdat%gauntarr, &
      !$acc        lmstart, cmt_nk, cmt_ikqpt)

         call timestart("gpu cpy in")
         !$acc wait
         call timestop("gpu cpy in")

         lm_0 = 0
         do iatm = 1,fi%atoms%nat
            itype = fi%atoms%itype(iatm)
            atom_phase = exp(-ImagUnit*tpi_const*dot_product(fi%kpts%bkf(:, iq), fi%atoms%taual(:, iatm)))
            atom_phase1 = sin(tpi_const*dot_product(fi%kpts%bkf(:, iq), fi%atoms%taual(:, iatm))) / sr2 ! This is not just the phase!!
            atom_phase2 = cos(tpi_const*dot_product(fi%kpts%bkf(:, iq), fi%atoms%taual(:, iatm))) / sr2 ! This is not just the phase!!

            iatm2 = fi%sym%invsatnr(iatm)
            IF(iatm2.EQ.0) iatm2 = iatm

            ioffset = sum((/((2*l + 1)*mpdata%num_radbasfn(l, itype), l=0, fi%hybinp%lcutm1(itype))/))


            IF(iatm2.LT.iatm) THEN ! iatm is the second of two atoms that are mapped onto each other by inversion symmetry

               DO l = 0, fi%hybinp%lcutm1(itype)
                  lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
               END DO
               CYCLE

            ELSE IF(iatm2.GT.iatm) THEN ! iatm is the first of two atoms that are mapped onto each other by inversion symmetry

               ! The default(shared) in the OMP part of the following loop is needed to avoid compilation issues on gfortran 7.5.
               DO l = 0, fi%hybinp%lcutm1(itype)
#ifdef _OPENACC
                  !$acc data copyin(l, iatm, iatm2, itype, lm_0, bandoi, bandof, psize, atom_phase, atom_phase1, atom_phase2, ioffset, ik, jsp)

                  !$acc parallel loop default(none) collapse(2)&
                  !$acc present(lmstart, cmt_ikqpt, cmt_nk, cprod,&
                  !$acc         l, iatm, iatm2, itype, lm_0, bandoi, bandof, psize, atom_phase, atom_phase1, atom_phase2, ioffset, ik, jsp, &
                  !$acc         mpdata, mpdata%num_radbasfn, mpdata%num_radfun_per_l, mpdata%l1, mpdata%l2, mpdata%n1, mpdata%n2,&
                  !$acc         hybdat, hybdat%prodm, hybdat%nbands, hybdat%nindxp1, hybdat%gauntarr)&
                  !$acc private(k,j,n,i,l1, l2, n1, n2, offdiag, lm1_0, lm2_0, lm, m, cscal, cscal2, add1, add2, ishift, lm1, m1, m2, lm2, lm1_cprod, lm2_cprod)
#else            
                  !$OMP PARALLEL DO default(shared) collapse(2) schedule(dynamic) & 
                  !$OMP private(k,j,n, n1, l1, n2, l2, offdiag, lm1_0, lm2_0, lm, m, cscal, cscal2, add1, add2, ishift, lm1, m1, m2, lm2, i, lm1_cprod, lm2_cprod)&
                  !$OMP shared(hybdat, bandoi, bandof, lmstart, lm_0, mpdata, cmt_ikqpt, cmt_nk, cprod, itype, l) &
                  !$OMP shared(iatm, iatm2, psize, atom_phase, atom_phase1, atom_phase2, ioffset, ik)
#endif
                  do k = 1, hybdat%nbands(ik,jsp) !This loop covers all bands that are to be used for the hybrid functionals calculation
                     do j = bandoi, bandof ! This loop only covers occupied bands
                        !$acc loop seq
                        DO n = 1, hybdat%nindxp1(l, itype) ! loop over basis-function products
                           ! don't call object funcktions in acc
                           l1 = mpdata%l1(n, l, itype) !
                           l2 = mpdata%l2(n, l, itype) ! current basis-function mpdatauct
                           n1 = mpdata%n1(n, l, itype) ! = bas(:,n1,l1,itype)*bas(:,n2,l2,itype) = b1*b2
                           n2 = mpdata%n2(n, l, itype) !

                           IF (mod(l1 + l2 + l, 2) == 0) THEN
                              offdiag = (l1 /= l2) .or. (n1 /= n2) ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                              !(leading to the same basis-function product)

                              lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                              lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                              lm = lm_0

                              !$acc loop seq
                              DO m = -l, l

                                 ishift = -2 * m * mpdata%num_radbasfn(l, itype)
                                 lm1_cprod = lm + (iatm - fi%atoms%firstAtom(itype))*ioffset
                                 lm2_cprod = lm + (iatm2 - fi%atoms%firstAtom(itype))*ioffset + ishift

                                 cscal = 0.0
                                 cscal2 = 0.0

                                 lm1 = lm1_0 + n1 ! go to lm index for m1=-l1
                                 !$acc loop seq
                                 DO m1 = -l1, l1
                                    m2 = m1 + m ! Gaunt condition -m1+m2-m=0
                                 
                                    IF (abs(m2) <= l2) THEN
                                       lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                       IF (abs(hybdat%gauntarr(1, l1, l2, l, m1, m)) > 1e-12) THEN
                                          cscal  = cscal  + hybdat%gauntarr(1, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm2, iatm))  * REAL(conjg(cmt_nk(k, lm1, iatm))) + &
                                                            hybdat%gauntarr(1, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm2, iatm2)) * REAL(conjg(cmt_nk(k, lm1, iatm2)))
                                          cscal2 = cscal2 + hybdat%gauntarr(1, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm2, iatm2)) * REAL(conjg(cmt_nk(k, lm1, iatm))) - &
                                                            hybdat%gauntarr(1, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm2, iatm))  * REAL(conjg(cmt_nk(k, lm1, iatm2)))
                                       END IF
                                    END IF

                                    m2 = m1 - m ! switch role of b1 and b2
                                    IF (abs(m2) <= l2 .and. offdiag) THEN
                                       lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                       IF (abs(hybdat%gauntarr(2, l1, l2, l, m1, m)) > 1e-12) THEN
                                          cscal  = cscal  + hybdat%gauntarr(2, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm1, iatm))  * REAL(conjg(cmt_nk(k, lm2, iatm))) + &
                                                            hybdat%gauntarr(2, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm1, iatm2)) * REAL(conjg(cmt_nk(k, lm2, iatm2)))
                                          cscal2 = cscal2 + hybdat%gauntarr(2, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm1, iatm2)) * REAL(conjg(cmt_nk(k, lm2, iatm))) - &
                                                            hybdat%gauntarr(2, l1, l2, l, m1, m) * REAL(cmt_ikqpt(j, lm1, iatm))  * REAL(conjg(cmt_nk(k, lm2, iatm2)))
                                       END IF
                                    END IF

                                    lm1 = lm1 + mpdata%num_radfun_per_l(l1, itype) ! go to lm start index for next m1-quantum number

                                 END DO  !m1

                                 add1 = cscal  * atom_phase2 + cscal2 * atom_phase1
                                 add2 = cscal2 * atom_phase2 - cscal  * atom_phase1

                                 !$acc loop seq
                                 DO i = 1, mpdata%num_radbasfn(l, itype)
                                    cprod(i + lm1_cprod, (j-bandoi+1) + (k-1)*psize) = cprod(i + lm1_cprod, (j-bandoi+1) + (k-1)*psize) + hybdat%prodm(i, n, l, itype) * add1
                                    cprod(i + lm2_cprod, (j-bandoi+1) + (k-1)*psize) = cprod(i + lm2_cprod, (j-bandoi+1) + (k-1)*psize) + hybdat%prodm(i, n, l, itype) * add2
                                 ENDDO
                                 lm = lm + mpdata%num_radbasfn(l, itype)
                              END DO ! m
                           ENDIF
                        END DO !n
                     enddo  !j
                  enddo !k
#ifdef _OPENACC
                  !$acc end parallel loop

                  !$acc end data 
#else
                  !$OMP END PARALLEL DO
#endif
                  lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
               END DO ! l loop

            ELSE ! (no inversion-symmetry mapping of atoms)

               ! The default(shared) in the OMP part of the following loop is needed to avoid compilation issues on gfortran 7.5.
               DO l = 0, fi%hybinp%lcutm1(itype)
#ifdef _OPENACC
                  !$acc data copyin(l, iatm, itype, lm_0, bandoi, bandof, psize, atom_phase, ik, jsp)

                  !$acc parallel loop default(none) collapse(2)&
                  !$acc present(lmstart, cmt_ikqpt, cmt_nk, cprod,&
                  !$acc         l, iatm, itype, lm_0, bandoi, bandof, psize, atom_phase, ik, jsp, &
                  !$acc         mpdata, mpdata%num_radbasfn, mpdata%num_radfun_per_l, mpdata%l1, mpdata%l2, mpdata%n1, mpdata%n2,&
                  !$acc         hybdat, hybdat%prodm, hybdat%nbands, hybdat%nindxp1, hybdat%gauntarr)&
                  !$acc private(k,j,n,i,l1, l2, n1, n2, offdiag, lm1_0, lm2_0, lm, m, cscal, lm1, m1, m2, lm2)
#else            
                  !$OMP PARALLEL DO default(shared) collapse(2) schedule(dynamic) & 
                  !$OMP private(k,j,n, n1, l1, n2, l2, offdiag, lm1_0, lm2_0, lm, m, cscal, lm1, m1, m2, lm2, i)&
                  !$OMP shared(hybdat, bandoi, bandof, lmstart, lm_0, mpdata, cmt_ikqpt, cmt_nk, cprod, itype, l) &
                  !$OMP shared(iatm, psize, atom_phase, ik)
#endif
                  do k = 1, hybdat%nbands(ik,jsp)
                     do j = bandoi, bandof 
                        !$acc loop seq
                        DO n = 1, hybdat%nindxp1(l, itype) ! loop over basis-function products
                           ! don't call object funcktions in acc
                           l1 = mpdata%l1(n, l, itype) !
                           l2 = mpdata%l2(n, l, itype) ! current basis-function mpdatauct
                           n1 = mpdata%n1(n, l, itype) ! = bas(:,n1,l1,itype)*bas(:,n2,l2,itype) = b1*b2
                           n2 = mpdata%n2(n, l, itype) !

                           IF (mod(l1 + l2 + l, 2) == 0) THEN
                              offdiag = (l1 /= l2) .or. (n1 /= n2) ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                              !(leading to the same basis-function product)

                              lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                              lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                              lm = lm_0
                              !$acc loop seq
                              DO m = -l, l
                                 cscal = 0.0

                                 lm1 = lm1_0 + n1 ! go to lm index for m1=-l1
                                 !$acc loop seq
                                 DO m1 = -l1, l1
                                    m2 = m1 + m ! Gaunt condition -m1+m2-m=0
                                 
                                    IF (abs(m2) <= l2) THEN
                                       lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                       IF (abs(hybdat%gauntarr(1, l1, l2, l, m1, m)) > 1e-12) THEN
                                          cscal = cscal + hybdat%gauntarr(1, l1, l2, l, m1, m) * cmt_ikqpt(j, lm2, iatm) * conjg(cmt_nk(k, lm1, iatm))
                                       END IF
                                    END IF

                                    m2 = m1 - m ! switch role of b1 and b2
                                    IF (abs(m2) <= l2 .and. offdiag) THEN
                                       lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                       IF (abs(hybdat%gauntarr(2, l1, l2, l, m1, m)) > 1e-12) THEN
                                          cscal = cscal + hybdat%gauntarr(2, l1, l2, l, m1, m) * cmt_ikqpt(j, lm1, iatm) * conjg(cmt_nk(k, lm2, iatm))
                                       END IF
                                    END IF

                                    lm1 = lm1 + mpdata%num_radfun_per_l(l1, itype) ! go to lm start index for next m1-quantum number

                                 END DO  !m1

                                 lm = lm_0 + (m + l)*mpdata%num_radbasfn(l, itype)
                                 !$acc loop seq
                                 DO i = 1, mpdata%num_radbasfn(l, itype)
                                    cprod(i + lm, (j-bandoi+1) + (k-1)*psize) = cprod(i + lm, (j-bandoi+1) + (k-1)*psize) + hybdat%prodm(i, n, l, itype)*cscal*atom_phase
                                 ENDDO
                              END DO
                           ENDIF
                        END DO !n
                     enddo  !j
                  enddo !k
#ifdef _OPENACC
                  !$acc end parallel loop

                  !$acc end data 
#else
                  !$OMP END PARALLEL DO
#endif
                  lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
               END DO ! l loop
            END IF
         END DO
      !$acc end data
      deallocate(cmt_ikqpt)

      call timestop("loop over l, l1, l2, n, n1, n2")
      call timestop("wavefproducts_noinv5 MT")
   end subroutine wavefproducts_noinv_MT
end module m_wavefproducts_noinv
