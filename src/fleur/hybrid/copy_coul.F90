module m_copy_coul
   use m_types
   use m_constants
   use m_glob_tofrom_loc
   USE m_types_mpimat
#ifdef CPP_MPI 
   use mpi 
#endif
   private 
   public :: copy_from_dense_to_sparse
contains
   subroutine copy_from_dense_to_sparse(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb(:)
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      call timestart("copy_from_dense_to_sparse")

      call copy_mt1_from_striped_to_sparse(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      call copy_mt2_from_striped_to_sparse(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      call copy_mt3_from_striped_to_sparse(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      call test_mt2_mt3(fi, fmpi, mpdata, ikpt, hybdat)
      call copy_residual_mt_contrib_atm(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      call copy_residual_mt_contrib_gpt(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      call copy_ir(fi, fmpi, mpdata, coulomb(ikpt), ikpt, hybdat)
      
      call timestop("copy_from_dense_to_sparse")
   end subroutine copy_from_dense_to_sparse

   subroutine copy_mt1_from_striped_to_sparse(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb(:)
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      real, allocatable    :: tmp_4r(:, :, :, :)
      complex, allocatable :: tmp_4c(:, :, :, :)
      integer :: indx1, sz, itype, ineq, l, i, ierr, ix_loc, n, pe_ix

      call timestart("copy_mt1")

      ! only one processor per k-point calculates MT convolution
      !
      ! store m-independent part of Coulomb matrix in MT spheres
      ! in coulomb_mt1(:mpdata%num_radbasfn(l,itype)-1,:mpdata%num_radbasfn(l,itype)-1,l,itype)
      !
      sz = maxval(mpdata%num_radbasfn) - 1
      if (fi%sym%invs) THEN
         allocate (tmp_4r(sz, sz, 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), source=0.0)
      else
         allocate (tmp_4c(sz, sz, 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), source=cmplx_0)
      end if
      indx1 = 0
      DO itype = 1, fi%atoms%ntype
         DO ineq = 1, fi%atoms%neq(itype)
            DO l = 0, fi%hybinp%lcutm1(itype)
               IF (ineq == 1) THEN
                  DO n = 1, mpdata%num_radbasfn(l, itype) - 1
                     do i = 1, mpdata%num_radbasfn(l, itype) - 1
                        call glob_to_loc(fmpi, indx1 + i, pe_ix, ix_loc)
                        if (fmpi%n_rank == pe_ix) then
                           if (fi%sym%invs) THEN
                              tmp_4r(n, i, l, itype) = real(coulomb(ikpt)%data_c(indx1 + n, ix_loc))
                           else
                              tmp_4c(n, i, l, itype) = real(coulomb(ikpt)%data_c(indx1 + n, ix_loc))
                           end if
                        end if
                     end do
                  END DO
               END IF

               indx1 = indx1 + (2*l + 1)*mpdata%num_radbasfn(l, itype)
            END DO
         END DO
      END do

      if (fi%sym%invs) THEN
#ifdef CPP_MPI
         call MPI_Reduce(tmp_4r, hybdat%coul(ikpt)%mt1_r, size(tmp_4r), MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, fmpi%sub_comm, ierr)
#else
         hybdat%coul(ikpt)%mt1_r = tmp_4r
#endif
         deallocate (tmp_4r)
      else
#ifdef CPP_MPI
         call MPI_Reduce(tmp_4c, hybdat%coul(ikpt)%mt1_c, size(tmp_4c), MPI_DOUBLE_COMPLEX, &
                         MPI_SUM, 0, fmpi%sub_comm, ierr)
#else
         hybdat%coul(ikpt)%mt1_c = tmp_4c
#endif
         deallocate (tmp_4c)
      end if
      call timestop("copy_mt1")
   end subroutine copy_mt1_from_striped_to_sparse


   subroutine copy_mt2_from_striped_to_sparse(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb(:)
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      integer :: indx1, itype, l, m, iatom, ic, ierr, info, ix, ix_loc, pe_ix, n
      real, allocatable    :: tmp_r(:,:,:,:)
      complex, allocatable :: tmp_c(:,:,:,:)

      call timestart("copy_mt2")
      if(fi%sym%invs) then
         allocate (tmp_r(maxval(mpdata%num_radbasfn) - 1, &
                        -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), &
                        0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info, source=0.0)
      else
         allocate (tmp_c(maxval(mpdata%num_radbasfn) - 1, &
                        -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), &
                        0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info, source=cmplx_0)
      endif
      if(info /= 0) call judft_error("can't alloc mt2_tmp")

      indx1 = 0
      do iatom = 1, fi%atoms%nat
         itype = fi%atoms%itype(iatom)
         DO l = 0, fi%hybinp%lcutm1(itype)
            DO M = -l, l
               ix = indx1 + mpdata%num_radbasfn(l, itype)
               call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
               if(pe_ix == fmpi%n_rank) then
                  if (fi%sym%invs) THEN
                     tmp_r(:mpdata%num_radbasfn(l, itype) - 1, M, l, iatom) &
                        = real(coulomb(ikpt)%data_c(indx1 + 1:indx1 + mpdata%num_radbasfn(l, itype) - 1, ix_loc))
                  else
                     tmp_c(:mpdata%num_radbasfn(l, itype) - 1, M, l, iatom) &
                        = coulomb(ikpt)%data_c(indx1 + 1:indx1 + mpdata%num_radbasfn(l, itype) - 1, ix_loc)
                  endif
               endif

               indx1 = indx1 + mpdata%num_radbasfn(l, itype)
            END DO
         END DO
      END DO

      ix = hybdat%nbasp + 1
      call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
      IF (ikpt == 1 .and. pe_ix == fmpi%n_rank) THEN
         !
         ! store the contribution of the G=0 plane wave with the MT l=0 functions in
         ! coulomb_mt2(:mpdata%num_radbasfn(l=0,itype),0,maxval(fi%hybinp%lcutm1)+1,iatom)
         !
         ic = 0
         do iatom = 1,fi%atoms%nat 
            itype = fi%atoms%itype(iatom)
            DO n = 1, mpdata%num_radbasfn(0, itype) - 1
               if (fi%sym%invs) THEN
                  tmp_r(n, 0, maxval(fi%hybinp%lcutm1) + 1, iatom) =  real(coulomb(ikpt)%data_c(ic + n, ix_loc))
               else
                  tmp_c(n, 0, maxval(fi%hybinp%lcutm1) + 1, iatom) = coulomb(ikpt)%data_c(ic + n, ix_loc)
               endif
            END DO
            ic = ic + SUM([((2*l + 1)*mpdata%num_radbasfn(l, itype), l=0, fi%hybinp%lcutm1(itype))])
         END DO
      endif 
      
      if (fi%sym%invs) THEN
#ifdef CPP_MPI
         call MPI_Reduce(tmp_r, hybdat%coul(ikpt)%mt2_r, size(tmp_r), MPI_DOUBLE_PRECISION, MPI_SUM, 0, fmpi%sub_comm, ierr)
#else
         hybdat%coul(ikpt)%mt2_r = tmp_r
#endif
         deallocate (tmp_r)
      else
#ifdef CPP_MPI
         call MPI_Reduce(tmp_c, hybdat%coul(ikpt)%mt2_c, size(tmp_c), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%sub_comm, ierr)
#else
         hybdat%coul(ikpt)%mt2_c = tmp_c
#endif
         deallocate (tmp_c)
      end if
      call timestop("copy_mt2")
   end subroutine copy_mt2_from_striped_to_sparse

   subroutine copy_mt3_from_striped_to_sparse(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      !
      ! store the contributions between the MT s-like functions at atom1 and
      ! and the constant function at a different atom2
      !
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb(:)
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      integer :: ic, iatom, itype, ishift, iatom1, ic1, ic2, itype1, ishift1, pe, loc, i, ierr, l, l1
      real, allocatable    :: tmp_r(:,:,:)
      complex, allocatable :: tmp_c(:,:,:)

      IF (ikpt == 1) THEN
         call timestart("copy_mt3")
         if(fi%sym%invs) then
            allocate(tmp_r(maxval(mpdata%num_radbasfn) - 1, fi%atoms%nat, fi%atoms%nat), source=0.0)
         else
            allocate(tmp_c(maxval(mpdata%num_radbasfn) - 1, fi%atoms%nat, fi%atoms%nat), source=cmplx_0)
         endif

         ic = 0
         do iatom = 1, fi%atoms%nat 
            itype = fi%atoms%itype(iatom)
            ishift = SUM([((2*l + 1)*mpdata%num_radbasfn(l, itype), l=0, fi%hybinp%lcutm1(itype))])
            ic1 = ic + mpdata%num_radbasfn(0, itype)

            ic2 = 0
            do iatom1 = 1,fi%atoms%nat
               itype1 = fi%atoms%itype(iatom1)
               ishift1 = SUM([((2*l1 + 1)*mpdata%num_radbasfn(l1, itype1), l1=0, fi%hybinp%lcutm1(itype1))])

               do i = 1,mpdata%num_radbasfn(0, itype1) - 1
                  call glob_to_loc(fmpi, ic2+i, pe, loc)
                  if(fmpi%n_rank == pe) then
                     IF (fi%sym%invs) THEN
                        tmp_r(i, iatom, iatom1) = real(coulomb(ikpt)%data_c(ic1, loc))
                     ELSE
                        tmp_c(i, iatom, iatom1) = CONJG(coulomb(ikpt)%data_c(ic1, loc))
                     ENDIF
                  endif
               enddo
               ic2 = ic2 + ishift1
            END DO
            ic = ic + ishift
         END DO

         if (fi%sym%invs) THEN
#ifdef CPP_MPI
            call MPI_Reduce(tmp_r, hybdat%coul(ikpt)%mt3_r, size(tmp_r), MPI_DOUBLE_PRECISION, MPI_SUM, 0, fmpi%sub_comm, ierr)
#else
            hybdat%coul(ikpt)%mt3_r = tmp_r
#endif
            deallocate (tmp_r)
         else
#ifdef CPP_MPI
            call MPI_Reduce(tmp_c, hybdat%coul(ikpt)%mt3_c, size(tmp_c), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%sub_comm, ierr)
#else
            hybdat%coul(ikpt)%mt3_c = tmp_c
#endif
            deallocate (tmp_c)
         end if
         call timestop("copy_mt3")
      endif ! ikpt == 1
   end subroutine copy_mt3_from_striped_to_sparse

   subroutine test_mt2_mt3(fi, fmpi, mpdata, ikpt, hybdat)
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      integer :: iatom, itype
      call timestart("test_mt2_mt3")
      if (fmpi%n_rank == 0 .and. ikpt == 1) then
         !test
         do iatom =1,fi%atoms%nat 
            itype = fi%atoms%itype(iatom)
            if (fi%sym%invs) THEN
               IF (MAXVAL(ABS(hybdat%coul(ikpt)%mt2_r(:mpdata%num_radbasfn(0, itype) - 1, 0, 0, iatom) &
                              - hybdat%coul(ikpt)%mt3_r(:mpdata%num_radbasfn(0, itype) - 1, iatom, iatom))) > 1E-08) &
                  call judft_error('coulombmatrix: coulomb_mt2 and coulomb_mt3 are inconsistent')

            else
               IF (MAXVAL(ABS(hybdat%coul(ikpt)%mt2_c(:mpdata%num_radbasfn(0, itype) - 1, 0, 0, iatom) &
                              - hybdat%coul(ikpt)%mt3_c(:mpdata%num_radbasfn(0, itype) - 1, iatom, iatom))) > 1E-08) &
                  call judft_error('coulombmatrix: coulomb_mt2 and coulomb_mt3 are inconsistent')
            end if
         END DO
      END IF
      call timestop("test_mt2_mt3")
   end subroutine test_mt2_mt3 

   subroutine copy_residual_mt_contrib_atm(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      !
      ! add the residual MT contributions, i.e. those functions with an moment,
      ! to the matrix coulomb_mtir, which is fully occupied
      !
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb(:)
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      integer :: igpt, indx1, indx2, indx3, indx4, itype, itype1, l, m, l1, m1
      integer :: iatom, iatom1, ierr, loc_4, pe_4, pe_ix, ix, ix_loc, ic, loc_from, i, tmp_idx
      complex :: tmp
      integer, allocatable :: loc_sizes(:), displs(:), loc_idx(:)
      complex, allocatable :: sendbuf(:), tmp_arr(:)
      
      call timestart("dbl iatom loop")
      ic = calc_ic(fi)

      allocate(loc_sizes(0:fmpi%n_size-1), displs(0:fmpi%n_size-1), loc_idx(0:fmpi%n_size-1))
      allocate(sendbuf(ic), tmp_arr(ic))
      indx1 = 0; indx2 = 0; indx3 = 0; indx4 = 0


      do iatom = 1, fi%atoms%nat 
         itype = fi%atoms%itype(iatom)
         DO l = 0, fi%hybinp%lcutm1(itype)
            DO M = -l, l
               indx1 = indx1 + 1
               indx3 = indx3 + mpdata%num_radbasfn(l, itype)


               loc_sizes = calc_loc_size_atom(fmpi, fi, mpdata, indx3)
               displs = calc_disp(loc_sizes)
               call assemble_sendbuf_atm(fi, fmpi, mpdata, coulomb, ikpt, indx3, sendbuf)
#ifdef CPP_MPI
               call MPI_Gatherv(sendbuf, loc_sizes(fmpi%n_rank), MPI_DOUBLE_COMPLEX, &
                               tmp_arr, loc_sizes, displs, MPI_DOUBLE_COMPLEX, 0, fmpi%sub_comm, ierr)
#else
               tmp_arr = sendbuf 
#endif

               if(fmpi%n_rank == 0) then
                  indx2 = 0
                  indx4 = 0
                  loc_idx = 0

                  do iatom1 = 1,fi%atoms%nat 
                     itype1 = fi%atoms%itype(iatom1)
                     DO l1 = 0, fi%hybinp%lcutm1(itype1)
                        DO m1 = -l1, l1
                           indx2 = indx2 + 1
                           indx4 = indx4 + mpdata%num_radbasfn(l1, itype1)
                           IF (indx4 >= indx3) then
                              call glob_to_loc(fmpi, indx4, pe_4, loc_4)
                              loc_idx(pe_4) = loc_idx(pe_4) + 1
                              IF (fi%sym%invs) THEN
                                 hybdat%coul(ikpt)%mtir%data_r(indx1, indx2) = real(tmp_arr(displs(pe_4) + loc_idx(pe_4)))
                              ELSE
                                 hybdat%coul(ikpt)%mtir%data_c(indx1, indx2) = tmp_arr(displs(pe_4) + loc_idx(pe_4))
                              ENDIF
                           endif
                        END DO
                     END DO
                  END DO
               endif !rank == 0
            enddo
         enddo
      enddo
      call timestop("dbl iatom loop")

   end subroutine copy_residual_mt_contrib_atm

   subroutine assemble_sendbuf_atm(fi, fmpi, mpdata, coulomb, ikpt, indx3, sendbuf)
      implicit none 
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb(:)
      integer, intent(in)               :: ikpt, indx3
      complex, intent(inout)            :: sendbuf(:)

      integer :: loc_idx, iatom1, itype1, l1, m1, indx4, pe_4, loc_4

      loc_idx = 0
      indx4 = 0

      do iatom1 = 1,fi%atoms%nat 
         itype1 = fi%atoms%itype(iatom1)
         DO l1 = 0, fi%hybinp%lcutm1(itype1)
            DO m1 = -l1, l1
               indx4 = indx4 + mpdata%num_radbasfn(l1, itype1)
               if (indx4 >= indx3) then
                  call glob_to_loc(fmpi, indx4, pe_4, loc_4)
                  if(pe_4 == fmpi%n_rank) then
                     loc_idx = loc_idx + 1
                     sendbuf(loc_idx) = coulomb(ikpt)%data_c(indx3, loc_4) 
                  endif
               endif
            enddo 
         enddo 
      enddo
   end subroutine assemble_sendbuf_atm

   subroutine copy_residual_mt_contrib_gpt(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb(:)
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      integer :: igpt, indx1, indx2, indx3, indx4, itype, itype1, l, m, l1, m1
      integer :: iatom, iatom1, ierr, loc_4, pe_4, pe_ix, ix, ix_loc, ic, loc_from, i, tmp_idx
      complex :: tmp
      complex, allocatable :: tmp_arr(:), sendbuf(:)
      integer, allocatable :: loc_sizes(:), displs(:), loc_froms(:)

      ic = calc_ic(fi)

      allocate(loc_sizes(0:fmpi%n_size-1), displs(0:fmpi%n_size-1), loc_froms(0:fmpi%n_size-1))
      loc_sizes = calc_loc_size_gpt(fmpi, hybdat, mpdata, ikpt)
      displs = calc_disp(loc_sizes)
      loc_froms = collect_loc_froms_gpt(fmpi, hybdat)
      allocate(sendbuf(loc_sizes(fmpi%n_rank)))

      allocate(tmp_arr(mpdata%n_g(ikpt)))
      call timestart("iatom igpt loop")
      indx1 = 0; indx3 = 0
      do iatom = 1, fi%atoms%nat 
         itype = fi%atoms%itype(iatom)
         DO l = 0, fi%hybinp%lcutm1(itype)
            DO M = -l, l
               indx1 = indx1 + 1
               indx3 = indx3 + mpdata%num_radbasfn(l, itype)

#ifdef CPP_MPI
               sendbuf = coulomb(ikpt)%data_c(indx3, loc_froms(fmpi%n_rank):)
               call MPI_Gatherv(sendbuf, loc_sizes(fmpi%n_rank), MPI_DOUBLE_COMPLEX, &
                               tmp_arr, loc_sizes, displs, MPI_DOUBLE_COMPLEX, 0, fmpi%sub_comm, ierr)
#else
               tmp_arr = coulomb(ikpt)%data_c(indx3, loc_froms(fmpi%n_rank):)
#endif

               if(fmpi%n_rank == 0) then
                  DO igpt = 1, mpdata%n_g(ikpt)
                     ix =  hybdat%nbasp + igpt
                     call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
                     tmp_idx = ix_loc - loc_froms(pe_ix) + 1 + displs(pe_ix)
                     IF (fi%sym%invs) THEN
                        hybdat%coul(ikpt)%mtir%data_r(indx1, ic + igpt) = real(tmp_arr(tmp_idx))
                     ELSE
                        hybdat%coul(ikpt)%mtir%data_c(indx1, ic + igpt) = tmp_arr(tmp_idx)
                     ENDIF
                  END DO
               endif

            END DO
         END DO
      END do
      call timestop("iatom igpt loop")

      call hybdat%coul(ikpt)%mtir%u2l()
      IF (indx1 /= ic) call judft_error('coulombmatrix: error index counting')
   end subroutine copy_residual_mt_contrib_gpt

   function calc_loc_size_gpt(fmpi, hybdat, mpdata, ikpt) result(loc_sizes)
      implicit none 
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      integer, intent(in)               :: ikpt
      integer :: loc_sizes(fmpi%n_size)
      integer :: loc_from, loc_to, my_size, ierr

      call range_from_glob_to_loc(fmpi, hybdat%nbasp + 1, loc_from)
      call range_to_glob_to_loc(fmpi, hybdat%nbasp + mpdata%n_g(ikpt), loc_to)

      my_size = loc_to - loc_from + 1
#ifdef CPP_MPI
      call MPI_Allgather(my_size, 1, MPI_INTEGER, loc_sizes, 1, MPI_INTEGER, fmpi%sub_comm, ierr)
#else 
      loc_sizes(1) = my_size 
#endif   
   end function calc_loc_size_gpt 

   function calc_loc_size_atom(fmpi, fi, mpdata, indx3) result(loc_size)
      implicit none 
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpi), intent(in)           :: fmpi
      type(t_mpdata), intent(in)        :: mpdata
      integer, intent(in)               :: indx3

      integer :: loc_size(0:fmpi%n_size-1)
      integer :: indx4, iatom1, itype1, l1, m1, loc_4, pe_4

      loc_size = 0
      indx4 = 0
      do iatom1 = 1,fi%atoms%nat 
         itype1 = fi%atoms%itype(iatom1)
         DO l1 = 0, fi%hybinp%lcutm1(itype1)
            DO m1 = -l1, l1
               indx4 = indx4 + mpdata%num_radbasfn(l1, itype1)
               IF (indx4 >= indx3) then
                  call glob_to_loc(fmpi, indx4, pe_4, loc_4)
                  loc_size(pe_4) = loc_size(pe_4) + 1
               endif
            enddo 
         enddo 
      enddo
   end function calc_loc_size_atom 

   function calc_disp(loc_sizes) result(displs)
      implicit NONE
      integer :: loc_sizes(:)
      integer :: displs(size(loc_sizes))
      integer :: i 

      displs = 0 
      do i = 2,size(loc_sizes)
         displs(i) = displs(i-1) + loc_sizes(i-1)
      end do
   end function calc_disp

   function collect_loc_froms_gpt(fmpi, hybdat) result(loc_froms)
      implicit none 
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      integer :: loc_froms(fmpi%n_size), ierr, loc_from

      call range_from_glob_to_loc(fmpi, hybdat%nbasp + 1, loc_from)

#ifdef CPP_MPI
      call MPI_Allgather(loc_from, 1, MPI_INTEGER, loc_froms, 1, MPI_INTEGER, fmpi%sub_comm, ierr)
#else 
      loc_froms(1) = loc_from 
#endif
   end function collect_loc_froms_gpt

   function calc_ic(fi) result(ic)
      implicit none 
      type(t_fleurinput), intent(in)    :: fi
      integer :: ic, iatom, itype, l

      ic = 0
      do iatom = 1,fi%atoms%nat 
         itype = fi%atoms%itype(iatom)
         ic = ic + (fi%hybinp%lcutm1(itype) + 1)**2
      END DO
   end function calc_ic

   subroutine copy_ir(fi, fmpi, mpdata, coulomb, ikpt, hybdat)
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      class(t_mat), intent(in)          :: coulomb
      integer, intent(in)               :: ikpt
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      integer :: ic, iatom, l, ix, iy, ix_loc, pe_ix, i, itype, ierr
      INTEGER:: blacs_desc(9), umap(1, 1), np
      real, allocatable    :: tmp(:)
      type(t_mat)          :: loc_cpy
      !
      ! add ir part to the matrix coulomb_mtir
      !

      ic = calc_ic(fi)
      call timestart("copy_ir")
      select type(coulomb)
      class is(t_mpimat)
#ifdef CPP_SCALAPACK
         if(fi%sym%invs) then
            call loc_cpy%alloc(.false., mpdata%n_g(ikpt), mpdata%n_g(ikpt))
            blacs_desc = [1, -1, loc_cpy%matsize1, loc_cpy%matsize2, loc_cpy%matsize1, loc_cpy%matsize2, 0, 0, loc_cpy%matsize1]
            umap(1, 1) = 0
            CALL BLACS_GET(coulomb%blacsdata%blacs_desc(2), 10, blacs_desc(2))
            CALL BLACS_GRIDMAP(blacs_desc(2), umap, 1, 1, 1)
            
            call pzgemr2d(mpdata%n_g(ikpt), mpdata%n_g(ikpt), &
                        coulomb%data_c, hybdat%nbasp + 1, hybdat%nbasp + 1, coulomb%blacsdata%blacs_desc, &
                        loc_cpy%data_c,1,  1,  blacs_desc, coulomb%blacsdata%blacs_desc(2))


            if(fmpi%n_rank == 0) then
               !$OMP parallel do default(shared) shared(mpdata, hybdat, loc_cpy, ic, ikpt) private(ix, iy) collapse(2)
               do ix = 1, mpdata%n_g(ikpt)
                  do iy = 1, mpdata%n_g(ikpt) 
                        hybdat%coul(ikpt)%mtir%data_r(ic + iy, ic + ix) = real(loc_cpy%data_c(iy, ix))
                  enddo 
               enddo
               !$OMP end parallel do
            endif
         else 
            blacs_desc = [1, -1, hybdat%coul(ikpt)%mtir%matsize1, hybdat%coul(ikpt)%mtir%matsize2, &
                        hybdat%coul(ikpt)%mtir%matsize1, hybdat%coul(ikpt)%mtir%matsize2, 0, 0, hybdat%coul(ikpt)%mtir%matsize1]
            umap(1, 1) = 0
            CALL BLACS_GET(coulomb%blacsdata%blacs_desc(2), 10, blacs_desc(2))
            CALL BLACS_GRIDMAP(blacs_desc(2), umap, 1, 1, 1)

            call pzgemr2d(mpdata%n_g(ikpt), mpdata%n_g(ikpt), &
                        coulomb%data_c, hybdat%nbasp + 1, hybdat%nbasp + 1, coulomb%blacsdata%blacs_desc, &
                        hybdat%coul(ikpt)%mtir%data_c, ic+1, ic+1,  blacs_desc, coulomb%blacsdata%blacs_desc(2))

         endif
#endif
      class is (t_mat)
         !$OMP parallel do default(shared) shared(mpdata, hybdat, loc_cpy, ic, fi, ikpt) private(ix, iy) collapse(2)
         do ix = 1, mpdata%n_g(ikpt)
            do iy = 1, mpdata%n_g(ikpt) 
               if(fi%sym%invs) then
                  hybdat%coul(ikpt)%mtir%data_r(ic + iy, ic + ix) = real(coulomb%data_c(hybdat%nbasp + iy, hybdat%nbasp + ix))
               else
                  hybdat%coul(ikpt)%mtir%data_c(ic + iy, ic + ix) = coulomb%data_c(hybdat%nbasp + iy, hybdat%nbasp + ix)
               endif
            enddo 
         enddo
         !$OMP end parallel do
      end select

      call timestop("copy_ir")
   end subroutine copy_ir
end module m_copy_coul
