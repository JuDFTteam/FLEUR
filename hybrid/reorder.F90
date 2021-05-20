MODULE m_reorder
#ifdef _OPENACC
   USE cublas
#define CPP_zswap cublasZswap
#define CPP_dswap cublasDswap
#else
#define CPP_zswap zswap
#define CPP_dswap dswap
#endif
   interface reorder
      module procedure reorder_real, reorder_cmplx
   end interface reorder
CONTAINS
   subroutine forw_order(atoms, lcutm, nindxm, new_order)
      USE m_types
      IMPLICIT NONE

      INTEGER, INTENT(IN)          :: lcutm(:), nindxm(0:, :)
      TYPE(t_atoms), INTENT(IN)    :: atoms
      integer, INTENT(INOUT)       :: new_order(:)
      
      INTEGER                  :: itype, ieq, indx1, indx2, l, n, m, info, i
      integer, allocatable     ::  tmp_order(:)

      new_order = [(i, i=1, size(new_order))]
      tmp_order = new_order

      indx1 = 0
      indx2 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, lcutm(itype)
               DO m = -l, l
                  DO n = 1, nindxm(l, itype) - 1
                     indx1 = indx1 + 1
                     indx2 = indx2 + 1
                     new_order(indx1) = tmp_order(indx2)
                  END DO
                  indx2 = indx2 + 1
               END DO
            END DO
         END DO
      END DO

      indx2 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, lcutm(itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + nindxm(l, itype)
                  new_order(indx1) = tmp_order(indx2)
               END DO
            END DO
         END DO
      END DO
   end subroutine forw_order

   subroutine back_order(atoms, lcutm, nindxm, new_order)
      use m_types 
      use m_judft
      implicit none 
      INTEGER, INTENT(IN)       :: lcutm(:), nindxm(0:, :)
      TYPE(t_atoms), INTENT(IN) :: atoms
      integer, INTENT(INOUT)    :: new_order(:)  
      
      INTEGER                   :: itype, ieq, indx1, indx2, l, n, m, info, i

      integer, allocatable     ::  tmp_order(:)

      new_order = [(i, i=1, size(new_order))]
      tmp_order = new_order

      indx1 = 0
      indx2 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, lcutm(itype)
               DO m = -l, l
                  DO n = 1, nindxm(l, itype) - 1
                     indx1 = indx1 + 1
                     indx2 = indx2 + 1
                     new_order(indx2) = tmp_order(indx1)
                  END DO
                  indx2 = indx2 + 1
               END DO
            END DO
         END DO
      END DO

      indx2 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, lcutm(itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + nindxm(l, itype)
                  new_order(indx2) = tmp_order(indx1)
               END DO
            END DO
         END DO
      END DO
   end subroutine back_order

   subroutine reorder_real(target_order, mat)
      implicit NONE 
      integer, intent(in)       :: target_order(:)
      REAL, INTENT(INOUT)       :: mat(:,:)

      integer, allocatable :: curr_order(:)
      integer :: i_tmp, i, j, sz_mat_1, sz_mat_2
      real    :: r_tmp

      sz_mat_1 = size(mat,1)
      sz_mat_2 = size(mat,2)
      curr_order = [(i, i=1, size(target_order))]

      do i = 1,size(mat, 1)
         if(curr_order(i) /= target_order(i)) then
            j = i + 1
            do while(target_order(i) /= curr_order(j))
               j = j + 1
            enddo

            i_tmp = curr_order(i)
            curr_order(i)   = curr_order(j)
            curr_order(j) = i_tmp

            !$acc host_data use_device(mat)
            call CPP_dswap(sz_mat_2, mat(i,1), sz_mat_1, mat(j,1), sz_mat_1)
            !$acc end host_data
         endif
      enddo
   end subroutine reorder_real

   subroutine reorder_cmplx(target_order, mat)
      implicit NONE 
      integer, intent(in)       :: target_order(:)
      complex, INTENT(INOUT)       :: mat(:,:)

      integer, allocatable :: curr_order(:)
      integer :: i_tmp, i, j, sz_mat_1, sz_mat_2
      complex    :: r_tmp

      sz_mat_1 = size(mat,1)
      sz_mat_2 = size(mat,2)
      curr_order = [(i, i=1, size(target_order))]

      do i = 1,size(mat,1)
         if(curr_order(i) /= target_order(i)) then
            j = i + 1
            do while(target_order(i) /= curr_order(j))
               j = j + 1
            enddo

            i_tmp = curr_order(i)
            curr_order(i)   = curr_order(j)
            curr_order(j) = i_tmp

            !$acc host_data use_device(mat)
            call CPP_zswap(size(mat,2), mat(i,1), sz_mat_1, mat(j,1), sz_mat_1)
            !$acc end host_data
         endif
      enddo
   end subroutine reorder_cmplx
END MODULE
