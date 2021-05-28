module m_spmm_noinv
   use iso_c_binding
   use m_spmm
#ifdef _OPENACC
   USE cublas
#define CPP_zgemm cublaszgemm
#define CPP_dgemm cublasdgemm
#define CPP_zgemv cublaszgemv
#define CPP_dgemv cublasdgemv
#define CPP_mtir_c mtir_tmp
#define CPP_mtir_r mtir_tmp
#else
#define CPP_zgemm zgemm
#define CPP_dgemm dgemm
#define CPP_zgemv zgemv
#define CPP_dgemv dgemv
#define CPP_mtir_c hybdat%coul(ikpt)%mtir%data_c
#define CPP_mtir_r hybdat%coul(ikpt)%mtir%data_r
#endif
contains
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
      complex, intent(inout)            :: mat_in(:,:)
      complex, intent(inout)            :: mat_out(:,:)

      integer :: n_vec, i_vec, ibasm, iatom, itype, ieq, l, m, n_size
      integer :: indx0, indx1, indx2, indx3, n, iatom1, ieq1, ishift, itype1
      integer :: ishift1, indx4, lm, idx1_start, idx3_start, ld_mt1_tmp
      integer :: iat2, it2, l2, iat, ierr, irank, i, sz_mtir, sz_in, sz_out, max_l_cut
      integer(C_SIZE_T) :: free_mem, tot_mem
      integer, allocatable :: new_order(:)
      complex, allocatable :: mt1_tmp(:,:,:,:), mt2_tmp(:,:,:,:), mt3_tmp(:,:,:), mat_in_line(:)
#ifdef _OPENACC
      complex, allocatable :: mtir_tmp(:,:)
#endif

      call timestart("spmm_noinvs")
      call timestart("copy mt2_c")
      mt2_tmp = hybdat%coul(ikpt)%mt2_c
      call timestop("copy mt2_c")

      sz_in  = size(mat_in, 1)
      sz_out  = size(mat_out, 1)
      n_vec = size(mat_in, 2)

      allocate(mat_in_line(size(mat_in,2)))

      call timestart("copyin gpu")
      !$acc data copyin(mt2_tmp) copy(mat_in) copyout(mat_out) create(mat_in_line)
         !$acc wait
         call timestop("copyin gpu")

         !$acc kernels present(mat_in_line, mat_in)
         mat_in_line = mat_in(hybdat%nbasp + 1, :)
         !$acc end kernels
         
         call timestart("reorder forw")
         allocate(new_order(size(mat_in,1)))
         call forw_order(fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, new_order)
         call reorder(new_order, mat_in)
         call timestop("reorder forw")

         ibasm = calc_ibasm(fi, mpdata)

         ! compute vecout for the indices from 0:ibasm
         call timestart("0 > ibasm: small matricies")
         call timestart("alloc&cpy mt1_tmp")
         allocate(mt1_tmp, mold=hybdat%coul(ikpt)%mt1_c, stat=ierr)
         ld_mt1_tmp = size(mt1_tmp,dim=2) ! special multiplication
         if(ierr /= 0) call judft_error("can't alloc mt1_tmp")
         call zcopy(size(mt1_tmp), hybdat%coul(ikpt)%mt1_c, 1, mt1_tmp, 1)
         call timestop("alloc&cpy mt1_tmp")


         !$acc kernels present(mat_out)
         mat_out = cmplx_0
         !$acc end kernels

         !$acc data copyin(mt1_tmp)
#ifndef _OPENACC
            !$OMP PARALLEL DO default(none) schedule(dynamic)&
            !$OMP private(iatom, itype, idx1_start, iat2, it2, l2, indx1, idx3_start, indx3)&
            !$OMP private(lm, l, m, n_size, i_vec)&
            !$OMP lastprivate(indx2)&
            !$OMP shared(ibasm, mat_in, hybdat, mat_out, fi, mpdata, n_vec, ikpt, ld_mt1_tmp, sz_out, sz_in, mt1_tmp, mt2_tmp)
#endif
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

                  !$acc host_data use_device(mt1_tmp, mat_in, mat_out)
                  call CPP_zgemm("N","N", n_size, n_vec, n_size, cmplx_1, mt1_tmp(1,1,l,itype), ld_mt1_tmp,&
                              mat_in(indx1,1), sz_in, cmplx_0, mat_out(indx1,1), sz_out)
                  !$acc end host_data

                  !$acc kernels present(mat_out, mt2_tmp, mat_in)
                  do i_vec = 1, n_vec
                     do i = 0, indx2-indx1
                        mat_out(indx1+i,i_vec) = mat_out(indx1+i,i_vec) + mt2_tmp(i+1, m, l, iatom) * mat_in(indx3, i_vec)
                     enddo
                  enddo
                  !$acc end kernels

                  indx1 = indx2
               END DO
            END DO
#ifndef _OPENACC
            !$OMP END PARALLEL DO
#endif
         !$acc end data
         !$acc wait
         deallocate(mt1_tmp)
         call timestop("0 > ibasm: small matricies")

         IF (indx2 /= ibasm) call judft_error('spmvec: error counting basis functions')

         IF (ikpt == 1) THEN
            call timestart("gamma point 1 noinv")
            call timestart("cpy mt3_tmp")
            allocate(mt3_tmp, mold=hybdat%coul(ikpt)%mt3_c, stat=ierr)
            if(ierr /= 0 ) call judft_error("can't alloc mt3_tmp")
            call zcopy(size(mt3_tmp), hybdat%coul(ikpt)%mt3_c, 1, mt3_tmp, 1)
            call timestop("cpy mt3_tmp")

            max_l_cut = maxval(fi%hybinp%lcutm1)
#ifdef _OPENACC
            !$acc data copyin(mt3_tmp)
#else
            !$OMP PARALLEL DO default(none) schedule(dynamic)&
            !$OMP private(iatom, itype, indx0, l, m, indx1, indx2, iatom1, indx3) &
            !$OMP private(indx4, i_vec, n_size, itype1, ishift1,ieq1) &
            !$OMP shared(fi, n_vec, mpdata, hybdat, ibasm, mat_out, mat_in, ikpt, mat_in_line, mt3_tmp, mt2_tmp, max_l_cut)
#endif
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
                        if (iatom /= iatom1) then
                           !$acc kernels present(mat_out, mt3_tmp, mat_in) default(none)
                           do i_vec = 1, n_vec
                              mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) + mt3_tmp(:n_size, iatom1, iatom)*mat_in(indx4, i_vec)
                           enddo
                           !$acc end kernels
                        endif
                     END DO
                     indx3 = indx3 + fi%atoms%neq(itype1)*ishift1
                  END DO
                  IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

                  n_size = mpdata%num_radbasfn(l, itype) - 1
                  !$acc kernels present(mat_out, mt2_tmp, mat_in_line) default(none)
                  do i_vec = 1, n_vec
                     mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) + mt2_tmp(:n_size, 0, max_l_cut + 1, iatom)*mat_in_line(i_vec)
                  enddo
                  !$acc end kernels
               END DO
#ifdef _OPENACC
            !$acc end data !(mt3_tmp)
#else
            !$OMP END PARALLEL DO
#endif
            call timestop("gamma point 1 noinv")
         END IF
         ! compute vecout for the index-range from ibasm+1:nbasm

         call timestart("calc indx1")
         indx1 = sum([(((2*l + 1)*fi%atoms%neq(itype), l=0, fi%hybinp%lcutm1(itype)), &
                     itype=1, fi%atoms%ntype)]) + mpdata%n_g(ikpt)
         call timestop("calc indx1")

#ifdef _OPENACC
         call timestart("copy mtir_tmp")
         allocate(mtir_tmp(hybdat%coul(ikpt)%mtir%matsize1, hybdat%coul(ikpt)%mtir%matsize2), stat=ierr)
         if(ierr /= 0) call judft_error("can't alloc mtir_tmp")
         call zlacpy("N", size(mtir_tmp,1), size(mtir_tmp,2), hybdat%coul(ikpt)%mtir%data_c, &
                     size(hybdat%coul(ikpt)%mtir%data_c,1), mtir_tmp, size(mtir_tmp,1))
         call timestop("copy mtir_tmp")
#endif

         call timestart("acc kernels")
         !$acc enter data copyin(mtir_tmp)
         if(conjg_mtir) then
            !$acc kernels present(mtir_tmp)
            CPP_mtir_c = conjg(CPP_mtir_c)
            !$acc end kernels
         endif
         call timestop("acc kernels")

         call timestart("ibasm+1->nbasm: zgemm")
         sz_mtir = size(CPP_mtir_c,1)

         !$acc host_data use_device(CPP_mtir_c, mat_in, mat_out)
         call CPP_zgemm("N", "N", indx1, n_vec, indx1, cmplx_1, CPP_mtir_c, sz_mtir, &
                     mat_in(ibasm + 1, 1), sz_in, cmplx_0, mat_out(ibasm + 1, 1), sz_out)
         !$acc end host_data
         !$acc exit data delete(CPP_mtir_c)
#ifdef _OPENACC
         deallocate(mtir_tmp)
#else       
         if(conjg_mtir) then
            CPP_mtir_c = conjg(CPP_mtir_c)
         endif
#endif
         !$acc wait
         call timestop("ibasm+1->nbasm: zgemm")

         call timestart("dot prod")
         !$acc kernels present(mt2_tmp)
         mt2_tmp = conjg(mt2_tmp)
         !$acc end kernels

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

                     !$acc host_data use_device(mat_in, mt2_tmp, mat_out)
                     call CPP_zgemv("T", n-1, n_vec, cmplx_1, mat_in(indx2,1), sz_in, mt2_tmp(1, m, l, iatom), 1, &
                     cmplx_1, mat_out(indx1,1), sz_out)
                     !$acc end host_data

                     indx2 = indx3
                  END DO

               END DO
            END DO
         END DO
         call timestop("dot prod")

         IF (ikpt == 1) THEN
            call timestart("gamma point 2 noinv")
            iatom = 0
            indx0 = 0

            max_l_cut = maxval(fi%hybinp%lcutm1)
            DO itype = 1, fi%atoms%ntype
               ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, fi%hybinp%lcutm1(itype))])
               DO ieq = 1, fi%atoms%neq(itype)
                  iatom = iatom + 1
                  indx1 = indx0 + 1
                  indx2 = indx1 + mpdata%num_radbasfn(0, itype) - 2
                  n_size = mpdata%num_radbasfn(0, itype) - 1

                  !$acc host_data use_device(mat_in, mt2_tmp, mat_out)
                  call CPP_zgemv("T", n_size, n_vec, cmplx_1, mat_in(indx1,1), sz_in, &
                     mt2_tmp(1,0,max_l_cut + 1, iatom), 1, cmplx_1, mat_out(hybdat%nbasp + 1, 1), sz_out)
                  !$acc end host_data
                  indx0 = indx0 + ishift
               END DO
            END DO

            !$acc data copyin(mt3_tmp)
               !$acc kernels present(mt3_tmp)
               mt3_tmp = conjg(mt3_tmp)
               !$acc end kernels
#ifndef _OPENACC
               !$OMP PARALLEL DO default(none) &
               !$OMP private(iatom, itype, indx1, iatom1, indx2, itype1, ishift1, indx3, indx4, n_size) &
               !$OMP shared(fi, mpdata, hybdat,mat_out, mat_in, ibasm, ikpt, n_vec, mt3_tmp, sz_out, sz_in)
#endif
               do iatom = 1, fi%atoms%nat 
                  itype = fi%atoms%itype(iatom)
                  indx1 = ibasm + sum([((fi%hybinp%lcutm1(fi%atoms%itype(iat)) + 1)**2, iat=1,iatom-1)]) + 1
                  iatom1 = 0
                  indx2 = 0
                  DO itype1 = 1, fi%atoms%ntype
                     ishift1 = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype1) - 1), l=0, fi%hybinp%lcutm1(itype1))])
                     DO ieq1 = 1, fi%atoms%neq(itype1)
                        iatom1 = iatom1 + 1
                        IF (iatom1 /= iatom) then
                           indx3 = indx2 + (ieq1 - 1)*ishift1 + 1
                           indx4 = indx3 + mpdata%num_radbasfn(0, itype1) - 2
                           n_size = mpdata%num_radbasfn(0, itype1) - 1

                           !$acc host_data use_device(mat_in, mt3_tmp, mat_out)
                           call CPP_zgemv("T", n_size, n_vec, cmplx_1, mat_in(indx3,1), sz_in, mt3_tmp(1, iatom, iatom1), 1, &
                                    cmplx_1, mat_out(indx1,1), sz_out)
                           !$acc end host_data
                        endif
                     END DO
                     indx2 = indx2 + fi%atoms%neq(itype1)*ishift1
                  END DO
               END DO
#ifndef _OPENACC
               !$OMP END PARALLEL DO
#endif
            !$acc end data !(mt3_tmp)
            deallocate(mt3_tmp) 
            call timestop("gamma point 2 noinv")
         END IF

         call timestart("reorder back")
         call back_order(fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, new_order)
         call reorder(new_order, mat_in)
         call reorder(new_order, mat_out)
         call timestop("reorder back")
      
         call timestart("copyout")
      !$acc end data !mt2_tmp, mat_in, mat_out
      !$acc wait
      call timestop("copyout")
      call timestop("spmm_noinvs")
   end subroutine spmm_noinvs
end module m_spmm_noinv