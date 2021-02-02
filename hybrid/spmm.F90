module m_spmm
   use iso_c_binding
#ifdef _OPENACC
   USE cublas
#define CPP_zgemm cublaszgemm
#define CPP_dgemm cublasdgemm
#define CPP_zgemv cublaszgemv
#define CPP_dgemv cublasdgemv
#else
#define CPP_zgemm zgemm
#define CPP_dgemm dgemm
#define CPP_zgemv zgemv
#define CPP_dgemv dgemv
#endif
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
      real, intent(inout)               :: mat_out(:,:)

      integer :: n_vec, i_vec, ibasm, iatom, itype, ieq, l, m, n_size, sz_mtir, sz_hlp, sz_out
      integer :: indx0, indx1, indx2, indx3, n, iatom1, ieq1, ishift, itype1
      integer :: ishift1, indx4, lm, iat2, it2, l2, idx1_start, idx3_start, iat, irank, ierr
      real, allocatable :: mat_hlp(:,:), mtir_tmp(:,:), mt2_tmp(:,:,:,:)

      call timestart("spmm_invs")
      allocate(mat_hlp(mat_in%matsize1, mat_in%matsize2), stat=ierr)
      if(ierr /= 0) call judft_error("can't alloc mat_hlp")
      mat_hlp = mat_in%data_r
      n_vec = mat_in%matsize2

      call timestart("reorder_forw")
      !$OMP PARALLEL DO default(none) &
      !$OMP private(i_vec) shared(n_vec, hybdat, ikpt, fi, mpdata, mat_hlp)
      do i_vec = 1, n_vec
         call reorder_forw(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_hlp(:, i_vec))
      enddo
      !$OMP end parallel do
      call timestop("reorder_forw")

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
            call dgemm("N","N", n_size, size(mat_hlp,2), n_size, 1.0, hybdat%coul(ikpt)%mt1_r(1,1,l,itype), size(hybdat%coul(ikpt)%mt1_r,dim=2),&
                        mat_hlp(indx1,1), size(mat_hlp,1), 0.0, mat_out(indx1,1), size(mat_out,1))

            do i_vec = 1, n_vec
               call daxpy(n_size, mat_hlp(indx3, i_vec), hybdat%coul(ikpt)%mt2_r(1,m,l,iatom), 1, mat_out(indx1,i_vec), 1)
            enddo

            indx1 = indx2
         END DO
      END DO
      !$OMP END PARALLEL DO
      call timestop("0 > ibasm: small matricies")

      IF (indx2 /= ibasm) call judft_error('spmm: error counting basis functions')

      IF (ikpt == 1) THEN
         call timestart("gamma point 1 inv")
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
                        mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt3_r(:n_size, iatom1, iatom)*mat_hlp(indx4, i_vec)
                     enddo
                  END DO
                  indx3 = indx3 + fi%atoms%neq(itype1)*ishift1
               END DO

               IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

               n_size = mpdata%num_radbasfn(l, itype) - 1
               do i_vec = 1, n_vec
                  mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt2_r(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom)*mat_in%data_r(indx3 + 1, i_vec)
               enddo
               indx0 = indx0 + ishift
            END DO
         END DO
         call timestop("gamma point 1 inv")
      END IF

      ! compute vecout for the index-range from ibasm+1:nbasm

      indx1 = sum([(((2*l + 1)*fi%atoms%neq(itype), l=0, fi%hybinp%lcutm1(itype)), &
                    itype=1, fi%atoms%ntype)]) + mpdata%n_g(ikpt)

      call timestart("ibasm+1 -> dgemm")

      call timestart("copy mtir_tmp")
      allocate(mtir_tmp(hybdat%coul(ikpt)%mtir%matsize1, hybdat%coul(ikpt)%mtir%matsize2), stat=ierr)
      if(ierr /= 0) call judft_error("can't alloc mtir_tmp")
      call dlacpy("N", size(mtir_tmp,1), size(mtir_tmp,2), hybdat%coul(ikpt)%mtir%data_r, &
                  size(hybdat%coul(ikpt)%mtir%data_r,1), mtir_tmp, size(mtir_tmp,1))
      call timestop("copy mtir_tmp")


      sz_mtir = size(mtir_tmp, 1)
      sz_hlp  = size(mat_hlp, 1)
      sz_out  = size(mat_out, 1)      
      
      !$acc enter data copyin(mtir_tmp, mat_hlp, mat_out)
      !$acc host_data use_device(mtir_tmp, mat_hlp, mat_out)  
      call CPP_dgemm("N", "N", indx1, n_vec, indx1, 1.0, mtir_tmp, sz_mtir, &
                 mat_hlp(ibasm + 1, 1), sz_hlp, 0.0, mat_out(ibasm + 1, 1), sz_out)
      !$acc end host_data
      !$acc exit data delete(mtir_tmp)
      deallocate(mtir_tmp)
      call timestop("ibasm+1 -> dgemm")

      call timestart("cpy mt2_tmp")
      mt2_tmp = hybdat%coul(ikpt)%mt2_r
      !$acc enter data copyin(mt2_tmp)
      call timestop("cpy mt2_tmp")

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

                  !$acc host_data use_device(mat_hlp, mt2_tmp, mat_out)
                  call CPP_dgemv("T", n-1, n_vec, 1.0, mat_hlp(indx2,1), sz_hlp, mt2_tmp(1, m, l, iatom), 1, &
                     1.0, mat_out(indx1,1), sz_out)
                  !$acc end host_data

                  indx2 = indx3
               END DO

            END DO
         END DO
      END DO
      call timestop("dot prod")

      !$acc exit data copyout(mat_out) delete(mt2_tmp, mat_hlp)


      IF (ikpt == 1) THEN
         call timestart("gamma point 2 inv")
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
                  mat_out(hybdat%nbasp + 1, i_vec) = mat_out(hybdat%nbasp + 1, i_vec) &
                                                            + dot_product(hybdat%coul(ikpt)%mt2_r(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom), mat_hlp(indx1:indx2, i_vec))
               enddo
               indx0 = indx0 + ishift
            END DO
         END DO

         !$OMP PARALLEL DO default(none) schedule(dynamic)&
         !$OMP private(iatom, itype, indx1, indx2, itype1, ishift1) &
         !$OMP private(ieq1, iatom1, indx3, indx4, n_size, i_vec) &
         !$OMP shared(fi, n_vec, mat_out, ibasm, mpdata, mat_hlp, hybdat, ikpt)
         do iatom = 1, fi%atoms%nat 
            itype = fi%atoms%itype(iatom)
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
                     mat_out(indx1, i_vec) = mat_out(indx1, i_vec) &
                                                      + dot_product(hybdat%coul(ikpt)%mt3_r(:n_size, iatom, iatom1), mat_hlp(indx3:indx4, i_vec))
                  enddo
               END DO
               indx2 = indx2 + fi%atoms%neq(itype1)*ishift1
            END DO
         END DO
         !$OMP END PARALLEL DO
         call timestop("gamma point 2 inv")

      END IF

      call timestart("reorder_back")
      !$OMP PARALLEL DO default(none) &
      !$OMP private(i_vec) shared(n_vec, hybdat, ikpt, fi, mpdata, mat_out)
      do i_vec = 1, n_vec
         call reorder_back(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_out(:, i_vec))
      enddo
      !$OMP END PARALLEL DO
      call timestop("reorder_back")
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
      type(t_hybdat), intent(inout)     :: hybdat
      integer, intent(in)               :: ikpt
      logical, intent(in)               :: conjg_mtir
      type(t_mat), intent(in)           :: mat_in
      complex, intent(inout)            :: mat_out(:,:)

      integer :: n_vec, i_vec, ibasm, iatom, itype, ieq, l, m, n_size
      integer :: indx0, indx1, indx2, indx3, n, iatom1, ieq1, ishift, itype1
      integer :: ishift1, indx4, lm, idx1_start, idx3_start
      integer :: iat2, it2, l2, iat, ierr, irank, i, sz_mtir, sz_hlp, sz_out
      integer(C_SIZE_T) :: free_mem, tot_mem
      complex, allocatable :: mat_hlp(:,:), mtir_tmp(:,:), mt2_tmp(:,:,:,:)

      call timestart("spmm_noinvs")
      allocate(mat_hlp(mat_in%matsize1, mat_in%matsize2), stat=ierr)
      if(ierr /= 0) call judft_error("can't alloc mat_hlp")
      mat_hlp = mat_in%data_c
      n_vec = mat_in%matsize2

      
      call timestart("reorder forw")
      !$OMP PARALLEL DO default(none) &
      !$OMP private(i_vec) shared(n_vec, hybdat, ikpt, fi, mpdata, mat_hlp)
      do i_vec = 1, n_vec
         call reorder_forw(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_hlp(:, i_vec))
      enddo
      !$OMP END PARALLEL DO
      call timestop("reorder forw")

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

            call zgemm("N","N", n_size, size(mat_hlp,2), n_size, cmplx_1, hybdat%coul(ikpt)%mt1_c(1,1,l,itype), size(hybdat%coul(ikpt)%mt1_c,dim=2),&
                        mat_hlp(indx1,1), size(mat_hlp,1), cmplx_0, mat_out(indx1,1), size(mat_out,1))

            do i_vec = 1, n_vec
               !call zaxpy(n_size, mat_hlp(indx3, i_vec), hybdat%coul(ikpt)%mt2_c(1,m,l,iatom), 1, mat_out(indx1,i_vec), 1)
               call zaxpy(n_size, mat_hlp(indx3, i_vec), hybdat%coul(ikpt)%mt2_c(1,m,l,iatom), 1, mat_out(indx1,i_vec), 1)
               !mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt2_c(:n_size, m, l, iatom)*mat_hlp(indx3, i_vec)
            enddo

            indx1 = indx2
         END DO
      END DO
      !$OMP END PARALLEL DO
      call timestop("0 > ibasm: small matricies")

      IF (indx2 /= ibasm) call judft_error('spmvec: error counting basis functions')

      IF (ikpt == 1) THEN
         call timestart("gamma point 1 noinv")
         !$OMP PARALLEL DO default(none) schedule(dynamic)&
         !$OMP private(iatom, itype, indx0, l, m, indx1, indx2, iatom1, indx3) &
         !$OMP private(indx4, i_vec, n_size, itype1, ishift1,ieq1) &
         !$OMP shared(fi, n_vec, mpdata, hybdat, ibasm, mat_out, mat_hlp, ikpt, mat_in)
         do iatom = 1,fi%atoms%nat
            itype = fi%atoms%itype(iatom)
            indx0 = 0
            do iat = 1,iatom-1
               indx0 = indx0 + sum([((2*l + 1)*(mpdata%num_radbasfn(l, fi%atoms%itype(iat)) - 1), l=0, fi%hybinp%lcutm1(fi%atoms%itype(iat)))])
            enddo
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
                     mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) + hybdat%coul(ikpt)%mt3_c(:n_size, iatom1, iatom)*mat_hlp(indx4, i_vec)
                  enddo
               END DO
               indx3 = indx3 + fi%atoms%neq(itype1)*ishift1
            END DO

            IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

            n_size = mpdata%num_radbasfn(l, itype) - 1
            do i_vec = 1, n_vec
               mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) &
                                                      + hybdat%coul(ikpt)%mt2_c(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom)*mat_in%data_c(indx3 + 1, i_vec)
            enddo
         END DO
         !$OMP END PARALLEL DO
         call timestop("gamma point 1 noinv")
      END IF
      ! compute vecout for the index-range from ibasm+1:nbasm

      call timestart("calc indx1")
      indx1 = sum([(((2*l + 1)*fi%atoms%neq(itype), l=0, fi%hybinp%lcutm1(itype)), &
                    itype=1, fi%atoms%ntype)]) + mpdata%n_g(ikpt)
      call timestop("calc indx1")


      call timestart("copy mtir_tmp")
      allocate(mtir_tmp(hybdat%coul(ikpt)%mtir%matsize1, hybdat%coul(ikpt)%mtir%matsize2), stat=ierr)
      if(ierr /= 0) call judft_error("can't alloc mtir_tmp")
      call zlacpy("N", size(mtir_tmp,1), size(mtir_tmp,2), hybdat%coul(ikpt)%mtir%data_c, &
                  size(hybdat%coul(ikpt)%mtir%data_c,1), mtir_tmp, size(mtir_tmp,1))
      call timestop("copy mtir_tmp")

      call timestart("acc kernels")
      !$acc enter data copyin(mtir_tmp, mat_hlp, mat_out)
      if(conjg_mtir) then
         !$acc kernels present(mtir_tmp)
         mtir_tmp = conjg(mtir_tmp)
         !$acc end kernels
      endif
      call timestop("acc kernels")

      call timestart("ibasm+1->nbasm: zgemm")
      sz_mtir = size(mtir_tmp,1)
      sz_hlp  = size(mat_hlp,1)
      sz_out  = size(mat_out, 1)

      !$acc host_data use_device(mtir_tmp, mat_hlp, mat_out)
      call CPP_zgemm("N", "N", indx1, n_vec, indx1, cmplx_1, mtir_tmp, sz_mtir, &
                    mat_hlp(ibasm + 1, 1), sz_hlp, cmplx_0, mat_out(ibasm + 1, 1), sz_out)
      !$acc end host_data
      !$acc exit data delete(mtir_tmp)
      deallocate(mtir_tmp)
      !$acc wait
      call timestop("ibasm+1->nbasm: zgemm")

      call timestart("copy mt2_c")
      mt2_tmp = hybdat%coul(ikpt)%mt2_c
      !$acc enter data copyin(mt2_tmp)
      !$acc kernels present(mt2_tmp)
      mt2_tmp = conjg(mt2_tmp)
      !$acc end kernels
      call timestop("copy mt2_c")

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

                  !$acc host_data use_device(mat_hlp, mt2_tmp, mat_out)
                  call CPP_zgemv("T", n-1, n_vec, cmplx_1, mat_hlp(indx2,1), sz_hlp, mt2_tmp(1, m, l, iatom), 1, &
                  cmplx_1, mat_out(indx1,1), sz_out)
                  !$acc end host_data

                  indx2 = indx3
               END DO

            END DO
         END DO
      END DO

      !$acc exit data copyout(mat_out) delete(mat_hlp, mt2_tmp)
      deallocate(mt2_tmp)
      call timestop("dot prod")

      IF (ikpt == 1) THEN
         call timestart("gamma point 2 noinv")
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
                  mat_out(hybdat%nbasp + 1, i_vec) = mat_out(hybdat%nbasp + 1, i_vec) &
                                                            + dot_product(hybdat%coul(ikpt)%mt2_c(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom), mat_hlp(indx1:indx2, i_vec))
               enddo
               indx0 = indx0 + ishift
            END DO
         END DO

         !$OMP PARALLEL DO default(none) &
         !$OMP private(iatom, itype, indx1, iatom1, indx2, itype1, ishift1, indx3, indx4, n_size) &
         !$OMP shared(fi, mpdata, hybdat,mat_out, mat_hlp, ibasm, ikpt, n_vec)
         do iatom = 1, fi%atoms%nat 
            itype = fi%atoms%itype(iatom)
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
                     mat_out(indx1, i_vec) = mat_out(indx1, i_vec) &
                                                      + dot_product(hybdat%coul(ikpt)%mt3_c(:n_size, iatom, iatom1), mat_hlp(indx3:indx4, i_vec))
                  enddo
               END DO
               indx2 = indx2 + fi%atoms%neq(itype1)*ishift1
            END DO
         END DO
         !$OMP END PARALLEL DO
         call timestop("gamma point 2 noinv")
      END IF

      call timestart("reorder back")
      !$OMP PARALLEL DO default(none) &
      !$OMP private(i_vec) shared(n_vec, hybdat, ikpt, fi, mpdata, mat_out)
      do i_vec = 1, n_vec
         call reorder_back(hybdat%nbasm(ikpt), fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, mat_out(:, i_vec))
      enddo
      !$OMP END PARALLEL DO
      call timestop("reorder back")
      call timestop("spmm_noinvs")
   end subroutine spmm_noinvs

   function calc_ibasm(fi, mpdata) result(ibasm)
      use m_types
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      integer :: ibasm, iatom, itype, ieq, l, m

      call timestart("calc_ibasm")
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
      call timestop("calc_ibasm")
   end function calc_ibasm
end module m_spmm
