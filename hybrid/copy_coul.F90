module m_copy_coul
   use m_types
   use m_constants
   use m_glob_tofrom_loc
#ifdef CPP_MPI 
   use mpi 
#endif
contains
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
end module m_copy_coul
