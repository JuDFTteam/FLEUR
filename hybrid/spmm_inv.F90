module m_spmm_inv
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
      real, intent(inout)               :: mat_in(:,:)
      real, intent(inout)               :: mat_out(:,:)

      integer :: n_vec, i_vec, ibasm, iatom, itype, ieq, l, m, n_size, sz_mtir, sz_hlp, sz_out, sz_mt1
      integer :: indx0, indx1, indx2, indx3, n, iatom1, ieq1, ishift, itype1, max_lcut_plus_1
      integer :: ishift1, indx4, lm, iat2, it2, l2, idx1_start, idx3_start, iat, irank, ierr
      real, allocatable :: mt1_tmp(:,:,:,:), mt2_tmp(:,:,:,:), mat_in_line(:), mt3_tmp(:,:,:)
      integer, allocatable :: new_order(:)
#ifdef _OPENACC 
      real, allocatable :: mtir_tmp(:,:)
#endif

      call timestart("spmm_invs")
      mat_in_line = mat_in(hybdat%nbasp + 1, :)
      
      n_vec = size(mat_in, 2)

      call timestart("reorder")
      allocate(new_order(size(mat_in,1)))
      call forw_order(fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, new_order)
      !$acc data copy(mat_in)
         call reorder(new_order, mat_in)
      !$acc end data
      call timestop("reorder")

      ibasm = calc_ibasm(fi, mpdata)


      call timestart("copies out of hydat")
      mt1_tmp = hybdat%coul(ikpt)%mt1_r
      mt2_tmp = hybdat%coul(ikpt)%mt2_r
      if(ikpt == 1 ) then
         mt3_tmp = hybdat%coul(ikpt)%mt3_r
      endif
      sz_mt1 = size(mt1_tmp,dim=2)
      sz_hlp  = size(mat_in, 1)
      sz_out  = size(mat_out, 1)     
      call timestop("copies out of hydat")


      !$acc data copyin(mat_in) copy(mat_out)
         !$acc data copyin(mt2_tmp)
            call timestart("0 > ibasm: small matricies")
            ! compute vecout for the indices from 0:ibasm
#ifndef _OPENACC
            !$OMP PARALLEL DO default(none) schedule(dynamic)&
            !$OMP private(iatom, itype, idx1_start, iat2, it2, l2, indx1, idx3_start, indx3)&
            !$OMP private(lm, l, m, n_size, i_vec)&
            !$OMP lastprivate(indx2)&
            !$OMP shared(ibasm, mat_in, hybdat, mat_out, fi, mpdata, n_vec, ikpt, mt2_tmp, sz_out, sz_hlp, sz_mt1, mt1_tmp)
#endif
            !$acc data copyin(mt1_tmp)
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
                     !$acc host_data use_device(mt1_tmp, mat_in, mat_out)  
                     call CPP_dgemm("N","N", n_size, n_vec, n_size, 1.0, mt1_tmp(1,1,l,itype), sz_mt1,&
                                 mat_in(indx1,1), sz_hlp, 0.0, mat_out(indx1,1), sz_out)
                     !$acc end host_data

                     !$acc kernels present(mat_out, mat_in, mt2_tmp)
                     do i_vec = 1, n_vec
                        mat_out(indx1:indx2,i_vec) = mat_out(indx1:indx2,i_vec) + mat_in(indx3, i_vec) * mt2_tmp(:n_size,m,l,iatom) 
                     enddo
                     !$acc end kernels
                     indx1 = indx2
                  END DO
               END DO
            !$acc end data
#ifndef _OPENACC
            !$OMP END PARALLEL DO
#endif
            call timestop("0 > ibasm: small matricies")

            IF (indx2 /= ibasm) call judft_error('spmm: error counting basis functions')

            IF (ikpt == 1) THEN
               !$acc data copyin(mt3_tmp, mat_in_line)
                  call timestart("gamma point 1 inv")
                  iatom = 0
                  indx0 = 0
#ifndef _OPENACC
                  !$OMP parallel do default(none) &
                  !$OMP private(iatom, itype, ishift, l, indx0, indx1, indx2, indx3, indx4, iatom1, itype1, ishift1, i_vec, n_size)&
                  !$OMP private(max_lcut_plus_1) shared(fi, mpdata, hybdat, mat_out, ibasm, n_vec, ikpt, mat_in, mat_in_line, mt2_tmp, mt3_tmp)
#endif
                  do iatom = 1, fi%atoms%nat 
                     itype = fi%atoms%itype(iatom)
                     ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, fi%hybinp%lcutm1(itype))])
                     l = 0

                     indx0 = 0
                     do iat = 1,iatom-1
                        indx0 = indx0 + sum([((2*l + 1)*(mpdata%num_radbasfn(l, fi%atoms%itype(iat)) - 1), l=0, fi%hybinp%lcutm1(fi%atoms%itype(iat)))])
                     enddo

                     indx1 = indx0 + 1
                     indx2 = indx1 + mpdata%num_radbasfn(l, itype) - 2

                     indx3 = ibasm
                     n_size = mpdata%num_radbasfn(l, itype) - 1
                     do iatom1 = 1,fi%atoms%nat 
                        itype1 = fi%atoms%itype(iatom1)
                        ishift1 = (fi%hybinp%lcutm1(itype1) + 1)**2
                        indx4 = indx3 + 1
                        IF (iatom /= iatom1) then
                           !$acc kernels present(mat_out, mat_in, mt3_tmp)
                           do i_vec = 1, n_vec
                              mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) &
                                 + mt3_tmp(:n_size, iatom1, iatom)*mat_in(indx4, i_vec)
                           enddo
                           !$acc end kernels
                        endif
                        indx3 = indx3 + ishift1
                     END DO

                     IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

                     n_size = mpdata%num_radbasfn(l, itype) - 1
                     max_lcut_plus_1 = maxval(fi%hybinp%lcutm1) + 1
                     !$acc kernels present(mat_out, mat_in_line, mt2_tmp)
                     do i_vec = 1, n_vec
                        mat_out(indx1:indx2, i_vec) = mat_out(indx1:indx2, i_vec) &
                           + mt2_tmp(:n_size, 0, max_lcut_plus_1, iatom)*mat_in_line(i_vec)
                     enddo
                     !$acc end kernels
                  END DO

#ifndef _OPENACC
                  !$OMP end parallel do
#endif
                  call timestop("gamma point 1 inv")
               !$acc end data !mt3_tmp
            END IF

            ! compute vecout for the index-range from ibasm+1:nbasm

            indx1 = sum([(((2*l + 1)*fi%atoms%neq(itype), l=0, fi%hybinp%lcutm1(itype)), &
                        itype=1, fi%atoms%ntype)]) + mpdata%n_g(ikpt)

            call timestart("ibasm+1 -> dgemm")
#ifdef _OPENACC
            call timestart("copy mtir_tmp")
            allocate(mtir_tmp(hybdat%coul(ikpt)%mtir%matsize1, hybdat%coul(ikpt)%mtir%matsize2), stat=ierr)
            if(ierr /= 0) call judft_error("can't alloc mtir_tmp")
            call dlacpy("N", size(mtir_tmp,1), size(mtir_tmp,2), hybdat%coul(ikpt)%mtir%data_r, &
                        size(hybdat%coul(ikpt)%mtir%data_r,1), mtir_tmp, size(mtir_tmp,1))
            call timestop("copy mtir_tmp")
#endif


            sz_mtir = size(CPP_mtir_r, 1)         
            !$acc data copyin(CPP_mtir_r)
               !$acc host_data use_device(CPP_mtir_r, mat_in, mat_out)  
               call CPP_dgemm("N", "N", indx1, n_vec, indx1, 1.0, CPP_mtir_r, sz_mtir, &
                        mat_in(ibasm + 1, 1), sz_hlp, 0.0, mat_out(ibasm + 1, 1), sz_out)
               !$acc end host_data
            !$acc end data ! CPP_mtir_r
#ifdef _OPENACC
            deallocate(mtir_tmp)
#endif
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

                        !$acc host_data use_device(mat_in, mt2_tmp, mat_out)
                        call CPP_dgemv("T", n-1, n_vec, 1.0, mat_in(indx2,1), sz_hlp, mt2_tmp(1, m, l, iatom), 1, &
                           1.0, mat_out(indx1,1), sz_out)
                        !$acc end host_data

                        indx2 = indx3
                     END DO

                  END DO
               END DO
            END DO
         !$acc end data ! mt2_tmp
         call timestop("dot prod")
      !$acc end data ! mat_in, mat_out



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
                                                   + dot_product(mt2_tmp(:n_size, 0, maxval(fi%hybinp%lcutm1) + 1, iatom), &
                                                                  mat_in(indx1:indx2, i_vec))
               enddo
               indx0 = indx0 + ishift
            END DO
         END DO

         !$OMP PARALLEL DO default(none) schedule(dynamic)&
         !$OMP private(iatom, itype, indx1, indx2, itype1, ishift1) &
         !$OMP private(ieq1, iatom1, indx3, indx4, n_size, i_vec) &
         !$OMP shared(fi, n_vec, mat_out, ibasm, mpdata, mat_in, hybdat, ikpt, mt3_tmp)
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
                                          + dot_product(mt3_tmp(:n_size, iatom, iatom1), &
                                                         mat_in(indx3:indx4, i_vec))
                  enddo
               END DO
               indx2 = indx2 + fi%atoms%neq(itype1)*ishift1
            END DO
         END DO
         !$OMP END PARALLEL DO
         call timestop("gamma point 2 inv")

      END IF

      call timestart("reorder") 
      call back_order(fi%atoms, fi%hybinp%lcutm1, mpdata%num_radbasfn, new_order)
      !$acc data copy(mat_in, mat_out)
         call reorder(new_order, mat_in)
         call reorder(new_order, mat_out)
      !$acc end data
      call timestop("reorder")
      call timestop("spmm_invs")
   end subroutine spmm_invs
end module m_spmm_inv