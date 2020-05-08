module m_spmm
! rewrite of spmvec to replace a sparse-matrix * vec multiplication by
! sparse-matrix * matrix
contains
   subroutine spmm_invs(fi, mpdata, hybdat, ikpt, mat_in, mat_out)
      use m_juDFT
      use m_types
      use m_reorder
      use m_calc_l_m_from_lm
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      type(t_hybdat), intent(in)        :: hybdat
      integer, intent(in)               :: ikpt
      type(t_mat), intent(in)           :: mat_in
      type(t_mat), intent(inout)        :: mat_out

      integer :: n_vec, i_vec, ibasm, iatom, itype, ieq, l, m, n_size
      integer :: indx0, indx1, indx2, indx3, n, iatom1, ieq1, ishift, itype1
      integer :: ishift1, indx4, lm, iat2, it2, l2, idx1_start, idx3_start, iat
      type(t_mat) :: mat_hlp, test_hlp, test_out

      call timestart("spmm_invs")
      call mat_hlp%init(mat_in)
      call mat_hlp%copy(mat_in, 1, 1)
      n_vec = mat_in%matsize2

      do i_vec = 1, n_vec
         call reorder_forw(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_hlp%data_r(:, i_vec))
      enddo

      ibasm = calc_ibasm(fi, mpdata)

      call timestart("0 > ibasm: small matricies")
      ! compute vecout for the indices from 0:ibasm
      !$OMP PARALLEL DO default(none) schedule(dynamic)&
      !$OMP private(iatom, itype, idx1_start, iat2, it2, l2, indx1, idx3_start, indx3)&
      !$OMP private(lm, l, m, n_size, i_vec)&
      !$OMP lastprivate(indx2)&
      !$OMP shared(ibasm, mat_hlp, hybdat, mat_out, fi, mpdata, n_vec, ikpt)
      do iatom = 1,fi%atoms%nat 
         itype = fi%atoms%itype(iatom)

         idx1_start = 0
         do iat2 =1,iatom-1
            it2 = fi%atoms%itype(iat2)
            do l2 = 0, fi%hybinp%lcutm1(it2)
               idx1_start = idx1_start + (mpdata%num_radbasfn(l2, it2)-1) * (2*l2+1)
            enddo
         enddo
         indx1 = idx1_start

         idx3_start = ibasm
         do iat2 = 1,iatom-1
            it2 = fi%atoms%itype(iat2)
            idx3_start = idx3_start + (fi%hybinp%lcutm1(it2)+1)**2
         enddo
         indx3 = idx3_start 

         do lm = 1,(fi%hybinp%lcutm1(itype)+1)**2
            call calc_l_m_from_lm(lm, l, m)
            indx1 = indx1 + 1
            indx2 = indx1 + mpdata%num_radbasfn(l, itype) - 2
            indx3 = indx3 + 1

            n_size = mpdata%num_radbasfn(l, itype) - 1
            call dgemm("N","N", n_size, mat_hlp%matsize2, n_size, 1.0, hybdat%coul(ikpt)%mt1_r(1,1,l,itype), size(hybdat%coul(ikpt)%mt1_r,dim=2),&
                        mat_hlp%data_r(indx1,1), mat_hlp%matsize1, 0.0, mat_out%data_r(indx1,1), mat_out%matsize1)

            do i_vec = 1, n_vec
               call daxpy(n_size, mat_hlp%data_r(indx3, i_vec), hybdat%coul(ikpt)%mt2_r(1,m,l,iatom), 1, mat_out%data_r(indx1,i_vec), 1)
            enddo

            indx1 = indx2
         END DO
      END DO
      !$OMP END PARALLEL DO
      call timestop("0 > ibasm: small matricies")

      IF (indx2 /= ibasm) call judft_error('spmm: error counting basis functions')

      IF (ikpt == 1) THEN
         iatom = 0
         indx0 = 0
         DO itype = 1, fi%atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, fi%hybinp%lcutm1(itype))])
            DO ieq = 1, fi%atoms%neq(itype)
               iatom = iatom + 1
               l = 0
               m = 0

               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(l, itype) - 2

               iatom1 = 0
               indx3 = ibasm
               n_size = mpdata%num_radbasfn(l, itype) - 1
               DO itype1 = 1, fi%atoms%ntype
                  ishift1 = (fi%hybinp%lcutm1(itype1) + 1)**2
                  DO ieq1 = 1, fi%atoms%neq(itype1)
                     iatom1 = iatom1 + 1
                     indx4 = indx3 + (ieq1 - 1)*ishift1 + 1
                     IF (iatom == iatom1) CYCLE
                     do i_vec = 1, n_vec
                        mat_out%data_r(indx1:indx2, i_vec) = mat_out%data_r(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt3_r(:n_size, iatom1, iatom)*mat_hlp%data_r(indx4, i_vec)
                     enddo
                  END DO
                  indx3 = indx3 + fi%atoms%neq(itype1)*ishift1
               END DO

               IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

               n_size = mpdata%num_radbasfn(l, itype) - 1
               do i_vec = 1, n_vec
                  mat_out%data_r(indx1:indx2, i_vec) = mat_out%data_r(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt2_r(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom)*mat_in%data_r(indx3 + 1, i_vec)
               enddo
               indx0 = indx0 + ishift
            END DO

         END DO
      END IF

      ! compute vecout for the index-range from ibasm+1:nbasm

      indx1 = sum([(((2*l + 1)*fi%atoms%neq(itype), l=0, fi%hybinp%lcutm1(itype)), &
                    itype=1, fi%atoms%ntype)]) + mpdata%n_g(ikpt)
      call timestart("ibasm+1 -> dgemm")
      ! call dgemm(transa, transb, m,   n,     k,  alpha,  a,          lda,   b,                         ldb,             beta,  c,                           ldc)
      call dgemm("N", "N", indx1, n_vec, indx1, 1.0, hybdat%coul(ikpt)%mtir_r, indx1, mat_hlp%data_r(ibasm + 1, 1), mat_hlp%matsize1, 0.0, mat_out%data_r(ibasm + 1, 1), mat_out%matsize1)
      call timestop("ibasm+1 -> dgemm")

      call timestart("dot prod")
      iatom = 0
      indx1 = ibasm; indx2 = 0; indx3 = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, fi%hybinp%lcutm1(itype)
               n = mpdata%num_radbasfn(l, itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + 1
                  indx3 = indx3 + n - 1

                  do i_vec = 1, n_vec
                     mat_out%data_r(indx1, i_vec) = mat_out%data_r(indx1, i_vec) + dot_product(hybdat%coul(ikpt)%mt2_r(:n - 1, m, l, iatom), mat_hlp%data_r(indx2:indx3, i_vec))
                  enddo
                  indx2 = indx3
               END DO

            END DO
         END DO
      END DO
      call timestop("dot prod")


      IF (ikpt == 1) THEN
         call timestart("gamma point 1 inv")
         iatom = 0
         indx0 = 0
         DO itype = 1, fi%atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, fi%hybinp%lcutm1(itype))])
            DO ieq = 1, fi%atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(0, itype) - 2
               n_size = mpdata%num_radbasfn(0, itype) - 1
               do i_vec = 1, n_vec
                  mat_out%data_r(hybdat%nbasp + 1, i_vec) = mat_out%data_r(hybdat%nbasp + 1, i_vec) &
                                                            + dot_product(hybdat%coul(ikpt)%mt2_r(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom), mat_hlp%data_r(indx1:indx2, i_vec))
               enddo
               indx0 = indx0 + ishift
            END DO
         END DO

         !$OMP PARALLEL DO default(none) schedule(dynamic)&
         !$OMP private(iatom, itype, ishift, indx1, indx2, itype1, ishift1) &
         !$OMP private(ieq1, iatom1, indx3, indx4, n_size, i_vec) &
         !$OMP shared(fi, n_vec, mat_out, ibasm, mpdata, mat_hlp, hybdat, ikpt)
         do iatom = 1, fi%atoms%nat 
            itype = fi%atoms%itype(iatom)
            ishift = (fi%hybinp%lcutm1(itype) + 1)**2
            indx1 = ibasm + sum([((fi%hybinp%lcutm1(fi%atoms%itype(iat)) + 1)**2, iat=1,iatom-1)]) + 1

            iatom1 = 0
            indx2 = 0
            DO itype1 = 1, fi%atoms%ntype
               ishift1 = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype1) - 1), l=0, fi%hybinp%lcutm1(itype1))])
               DO ieq1 = 1, fi%atoms%neq(itype1)
                  iatom1 = iatom1 + 1
                  IF (iatom1 == iatom) CYCLE

                  indx3 = indx2 + (ieq1 - 1)*ishift1 + 1
                  indx4 = indx3 + mpdata%num_radbasfn(0, itype1) - 2

                  n_size = mpdata%num_radbasfn(0, itype1) - 1
                  do i_vec = 1, n_vec
                     mat_out%data_r(indx1, i_vec) = mat_out%data_r(indx1, i_vec) &
                                                      + dot_product(hybdat%coul(ikpt)%mt3_r(:n_size, iatom, iatom1), mat_hlp%data_r(indx3:indx4, i_vec))
                  enddo
               END DO
               indx2 = indx2 + fi%atoms%neq(itype1)*ishift1
            END DO
         END DO
         !$OMP END PARALLEL DO

         call timestop("gamma point 1 inv")
      END IF

      do i_vec = 1, n_vec
         call reorder_back(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_out%data_r(:, i_vec))
      enddo
      call timestop("spmm_invs")
   end subroutine spmm_invs

   subroutine spmm_noinvs(fi, mpdata, hybdat, ikpt, conjg_mtir, mat_in, mat_out)
      use m_juDFT
      use m_types
      use m_reorder
      use m_constants
      use m_calc_l_m_from_lm
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      type(t_hybdat), intent(in)        :: hybdat
      integer, intent(in)               :: ikpt
      logical, intent(in)               :: conjg_mtir
      type(t_mat), intent(in)           :: mat_in
      type(t_mat), intent(inout)        :: mat_out

      integer :: n_vec, i_vec, ibasm, iatom, itype, ieq, l, m, n_size
      integer :: indx0, indx1, indx2, indx3, n, iatom1, ieq1, ishift, itype1
      integer :: ishift1, indx4, iatom2, l1, lm, idx1_start, idx3_start
      integer :: iat2, it2, l2
      type(t_mat) :: mat_hlp, test_hlp, test_out

      call timestart("spmm_noinvs")
      call mat_hlp%init(mat_in)
      call mat_hlp%copy(mat_in, 1, 1)
      n_vec = mat_in%matsize2

      do i_vec = 1, n_vec
         call reorder_forw(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_hlp%data_c(:, i_vec))
      enddo

      ibasm = calc_ibasm(fi, mpdata)

      ! compute vecout for the indices from 0:ibasm

      call timestart("0 > ibasm: small matricies")
      !$OMP PARALLEL DO default(none) schedule(dynamic)&
      !$OMP private(iatom, itype, idx1_start, iat2, it2, l2, indx1, idx3_start, indx3)&
      !$OMP private(lm, l, m, n_size, i_vec)&
      !$OMP lastprivate(indx2)&
      !$OMP shared(ibasm, mat_hlp, hybdat, mat_out, fi, mpdata, n_vec, ikpt)
      do iatom = 1, fi%atoms%nat
         itype = fi%atoms%itype(iatom)

         idx1_start = 0
         do iat2 =1,iatom-1
            it2 = fi%atoms%itype(iat2)
            do l2 = 0, fi%hybinp%lcutm1(it2)
               idx1_start = idx1_start + (mpdata%num_radbasfn(l2, it2)-1) * (2*l2+1)
            enddo
         enddo
         indx1 = idx1_start

         idx3_start = ibasm
         do iat2 = 1,iatom-1
            it2 = fi%atoms%itype(iat2)
            idx3_start = idx3_start + (fi%hybinp%lcutm1(it2)+1)**2
         enddo
         indx3 = idx3_start 
         do lm = 1, (fi%hybinp%lcutm1(itype) + 1)**2
            call calc_l_m_from_lm(lm, l, m)
            indx1 = indx1 + 1
            indx2 = indx1 + mpdata%num_radbasfn(l, itype) - 2
            indx3 = indx3 + 1

            n_size = mpdata%num_radbasfn(l, itype) - 1

            call zgemm("N","N", n_size, mat_hlp%matsize2, n_size, cmplx_1, hybdat%coul(ikpt)%mt1_c(1,1,l,itype), size(hybdat%coul(ikpt)%mt1_c,dim=2),&
                        mat_hlp%data_c(indx1,1), mat_hlp%matsize1, cmplx_0, mat_out%data_c(indx1,1), mat_out%matsize1)

            do i_vec = 1, n_vec
               !call zaxpy(n_size, mat_hlp%data_c(indx3, i_vec), hybdat%coul(ikpt)%mt2_c(1,m,l,iatom), 1, mat_out%data_c(indx1,i_vec), 1)
               call zaxpy(n_size, mat_hlp%data_c(indx3, i_vec), hybdat%coul(ikpt)%mt2_c(1,m,l,iatom), 1, mat_out%data_c(indx1,i_vec), 1)
               !mat_out%data_c(indx1:indx2, i_vec) = mat_out%data_c(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt2_c(:n_size, m, l, iatom)*mat_hlp%data_c(indx3, i_vec)
            enddo

            indx1 = indx2
         END DO
      END DO
      !$OMP END PARALLEL DO
      call timestop("0 > ibasm: small matricies")

      IF (indx2 /= ibasm) call judft_error('spmvec: error counting basis functions')

      IF (ikpt == 1) THEN
         call timestart("gamma point 1")
         iatom = 0
         indx0 = 0
         DO itype = 1, fi%atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, fi%hybinp%lcutm1(itype))])
            DO ieq = 1, fi%atoms%neq(itype)
               iatom = iatom + 1
               l = 0
               m = 0

               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(l, itype) - 2

               iatom1 = 0
               indx3 = ibasm
               n_size = mpdata%num_radbasfn(l, itype) - 1
               DO itype1 = 1, fi%atoms%ntype
                  ishift1 = (fi%hybinp%lcutm1(itype1) + 1)**2
                  DO ieq1 = 1, fi%atoms%neq(itype1)
                     iatom1 = iatom1 + 1
                     indx4 = indx3 + (ieq1 - 1)*ishift1 + 1
                     IF (iatom == iatom1) CYCLE
                     do i_vec = 1, n_vec
                        mat_out%data_c(indx1:indx2, i_vec) = mat_out%data_c(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt3_c(:n_size, iatom1, iatom)*mat_hlp%data_c(indx4, i_vec)
                     enddo
                  END DO
                  indx3 = indx3 + fi%atoms%neq(itype1)*ishift1
               END DO

               IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

               n_size = mpdata%num_radbasfn(l, itype) - 1
               do i_vec = 1, n_vec
                  mat_out%data_c(indx1:indx2, i_vec) = mat_out%data_c(indx1:indx2, i_vec) &
                                                       + hybdat%coul(ikpt)%mt2_c(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom)*mat_in%data_c(indx3 + 1, i_vec)
               enddo

               indx0 = indx0 + ishift
            END DO

         END DO
         call timestop("gamma point 1")
      END IF
      ! compute vecout for the index-range from ibasm+1:nbasm

      indx1 = sum([(((2*l + 1)*fi%atoms%neq(itype), l=0, fi%hybinp%lcutm1(itype)), &
                    itype=1, fi%atoms%ntype)]) + mpdata%n_g(ikpt)
      call timestart("ibasm+1->nbasm: zgemm")
      if(conjg_mtir) then
         call zgemm("N", "N", indx1, n_vec, indx1, cmplx_1, conjg(hybdat%coul(ikpt)%mtir_c), indx1, &
                    mat_hlp%data_c(ibasm + 1, 1), mat_hlp%matsize1, cmplx_0, mat_out%data_c(ibasm + 1, 1), mat_out%matsize1)
      else
         call zgemm("N", "N", indx1, n_vec, indx1, cmplx_1, hybdat%coul(ikpt)%mtir_c, indx1, &
                    mat_hlp%data_c(ibasm + 1, 1), mat_hlp%matsize1, cmplx_0, mat_out%data_c(ibasm + 1, 1), mat_out%matsize1)
      endif
      call timestop("ibasm+1->nbasm: zgemm")

      call timestart("dot prod")
      iatom = 0
      indx1 = ibasm; indx2 = 0; indx3 = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, fi%hybinp%lcutm1(itype)
               n = mpdata%num_radbasfn(l, itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + 1
                  indx3 = indx3 + n - 1

                  do i_vec = 1, n_vec
                     mat_out%data_c(indx1, i_vec) = mat_out%data_c(indx1, i_vec) + dot_product(hybdat%coul(ikpt)%mt2_c(:n - 1, m, l, iatom), mat_hlp%data_c(indx2:indx3, i_vec))
                  enddo
                  indx2 = indx3
               END DO

            END DO
         END DO
      END DO
      call timestop("dot prod")

      IF (ikpt == 1) THEN
         call timestart("gamma point 2")
         iatom = 0
         indx0 = 0
         DO itype = 1, fi%atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, fi%hybinp%lcutm1(itype))])
            DO ieq = 1, fi%atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(0, itype) - 2
               n_size = mpdata%num_radbasfn(0, itype) - 1
               do i_vec = 1, n_vec
                  mat_out%data_c(hybdat%nbasp + 1, i_vec) = mat_out%data_c(hybdat%nbasp + 1, i_vec) &
                                                            + dot_product(hybdat%coul(ikpt)%mt2_c(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom), mat_hlp%data_c(indx1:indx2, i_vec))
               enddo
               indx0 = indx0 + ishift
            END DO
         END DO

         iatom = 0
         indx0 = ibasm
         DO itype = 1, fi%atoms%ntype
            ishift = (fi%hybinp%lcutm1(itype) + 1)**2
            DO ieq = 1, fi%atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1

               iatom1 = 0
               indx2 = 0
               DO itype1 = 1, fi%atoms%ntype
                  ishift1 = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype1) - 1), l=0, fi%hybinp%lcutm1(itype1))])
                  DO ieq1 = 1, fi%atoms%neq(itype1)
                     iatom1 = iatom1 + 1
                     IF (iatom1 == iatom) CYCLE

                     indx3 = indx2 + (ieq1 - 1)*ishift1 + 1
                     indx4 = indx3 + mpdata%num_radbasfn(0, itype1) - 2
                     n_size = mpdata%num_radbasfn(0, itype1) - 1
                     do i_vec = 1, n_vec
                        mat_out%data_c(indx1, i_vec) = mat_out%data_c(indx1, i_vec) &
                                                       + dot_product(hybdat%coul(ikpt)%mt3_c(:n_size, iatom, iatom1), mat_hlp%data_c(indx3:indx4, i_vec))
                     enddo
                  END DO
                  indx2 = indx2 + fi%atoms%neq(itype1)*ishift1
               END DO
               indx0 = indx0 + ishift
            END DO
         END DO
         IF (indx0 /= hybdat%nbasp) call judft_error('spmvec: error index counting (indx0)')
         call timestop("gamma point 2")
      END IF

      do i_vec = 1, n_vec
         call reorder_back(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_out%data_c(:, i_vec))
      enddo
      call timestop("spmm_noinvs")
   end subroutine spmm_noinvs

   subroutine spmv_wrapper_inv(fi, mpdata, hybdat, ikpt, vecin, vecout)
      use m_types
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      type(t_hybdat), intent(in)        :: hybdat
      integer, intent(in)               :: ikpt
      REAL, INTENT(IN) ::  vecin(:)!(hybdat%nbasm)
      REAL, INTENT(INOUT)::  vecout(:)!(hybdat%nbasm)

      type(t_mat) :: mat_in, mat_out

      call mat_in%alloc(.true., size(vecin), 1)
      call mat_out%alloc(.true., size(vecout), 1)

      mat_in%data_r(:, 1) = vecin
      mat_out%data_r(:, 1) = vecout
      call spmm_invs(fi, mpdata, hybdat, ikpt, mat_in, mat_out)
      vecout = mat_out%data_r(:, 1)
   end subroutine spmv_wrapper_inv

   subroutine spmv_wrapper_noinv(fi, mpdata, hybdat, ikpt, conjg_mtir, vecin, vecout)
      use m_types
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      type(t_hybdat), intent(in)        :: hybdat
      integer, intent(in)               :: ikpt
      logical, intent(in)               :: conjg_mtir
      complex, INTENT(IN)               :: vecin(:)!(hybdat%nbasm)
      complex, INTENT(INOUT)            :: vecout(:)!(hybdat%nbasm)

      type(t_mat) :: mat_in, mat_out

      call mat_in%alloc(.false., size(vecin), 1)
      call mat_out%alloc(.false., size(vecout), 1)

      mat_in%data_c(:, 1) = vecin
      mat_out%data_c(:, 1) = vecout
      call spmm_noinvs(fi, mpdata, hybdat, ikpt, conjg_mtir, mat_in, mat_out)
      vecout = mat_out%data_c(:, 1)
   end subroutine spmv_wrapper_noinv

   function calc_ibasm(fi, mpdata) result(ibasm)
      use m_types
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      integer :: ibasm, iatom, itype, ieq, l, m

      ibasm = 0
      iatom = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, fi%hybinp%lcutm1(itype)
               DO m = -l, l
                  ibasm = ibasm + mpdata%num_radbasfn(l, itype) - 1
               END DO
            END DO
         END DO
      END DO
   end function calc_ibasm
end module m_spmm
