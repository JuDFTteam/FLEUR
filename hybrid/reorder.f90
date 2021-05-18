MODULE m_reorder
   interface reorder_forw
      module procedure reorder_forw_real, reorder_forw_cmplx
   end interface reorder_forw

   interface reorder_back
      module procedure reorder_back_real, reorder_back_cmplx
   end interface reorder_back

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

   ! subroutine back_order(atoms, lcutm, nindxm, new_order)
   !    use m_types 
   !    use m_judft
   !    implicit none 
   !    INTEGER, INTENT(IN)       :: lcutm(:), nindxm(0:, :)
   !    TYPE(t_atoms), INTENT(IN) :: atoms
   !    REAL, INTENT(INOUT)       :: vec_r(nbasm)
      
   !    INTEGER                   :: itype, ieq, indx1, indx2, l, n, m, info
   !    integer, INTENT(INOUT)    :: new_order(:)  
   ! end subroutine back_order

   subroutine reorder_forw_real(target_order, mat)
      implicit NONE 
      integer, intent(in)       :: target_order(:)
      REAL, INTENT(INOUT)       :: mat(:,:)

      integer, allocatable :: curr_order(:)
      integer :: i_tmp, i, j, ld_mat 
      real    :: r_tmp

      ld_mat = size(mat,1)
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

            call dswap(size(mat,2), mat(i,1), ld_mat, mat(j,1), ld_mat)
         endif
      enddo
   end subroutine reorder_forw_real

   subroutine reorder_forw_cmplx(target_order, mat)
      implicit NONE 
      integer, intent(in)       :: target_order(:)
      complex, INTENT(INOUT)       :: mat(:,:)

      integer, allocatable :: curr_order(:)
      integer :: i_tmp, i, j, ld_mat
      complex    :: r_tmp

      ld_mat = size(mat,1)
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

            call zswap(size(mat,2), mat(i,1), ld_mat, mat(j,1), ld_mat)
         endif
      enddo
   end subroutine reorder_forw_cmplx

   subroutine reorder_back_real(nbasm, atoms, lcutm, nindxm, vec_r)
      use m_types 
      use m_judft
      implicit none 
      INTEGER, INTENT(IN)       :: nbasm, lcutm(:), nindxm(0:, :)
      TYPE(t_atoms), INTENT(IN) :: atoms
      REAL, INTENT(INOUT)       :: vec_r(nbasm)
      
      INTEGER               :: itype, ieq, indx1, indx2, l, n, m, info
      REAL, allocatable     ::  vechlp_r(:)      
      
      allocate(vechlp_r(nbasm), source=0.0, stat=info)
      if(info /= 0) call judft_error("can't allocate vechlp_r")
      vechlp_r = vec_r

      indx1 = 0
      indx2 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, lcutm(itype)
               DO m = -l, l
                  DO n = 1, nindxm(l, itype) - 1
                     indx1 = indx1 + 1
                     indx2 = indx2 + 1
                     vec_r(indx2) = vechlp_r(indx1)
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
                  vec_r(indx2) = vechlp_r(indx1)
               END DO
            END DO
         END DO
      END DO
   end subroutine reorder_back_real

   subroutine reorder_back_cmplx(nbasm, atoms, lcutm, nindxm, vec_c)
      USE m_types
      USE m_juDFT
      use m_constants, only: cmplx_0
      IMPLICIT NONE

      INTEGER, INTENT(IN)       :: nbasm, lcutm(:), nindxm(0:, :)
      TYPE(t_atoms), INTENT(IN) :: atoms
      complex, INTENT(INOUT)    :: vec_c(nbasm)
      
      INTEGER               :: itype, ieq, indx1, indx2, l, n, m, info
      complex, allocatable  :: vechlp_c(:)

      allocate(vechlp_c(nbasm), source=cmplx_0, stat=info)
      if(info /= 0) call judft_error("can't allocate vechlp_c")
      vechlp_c = vec_c

      indx1 = 0
      indx2 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, lcutm(itype)
               DO m = -l, l
                  DO n = 1, nindxm(l, itype) - 1
                     indx1 = indx1 + 1
                     indx2 = indx2 + 1
                     vec_c(indx2) = vechlp_c(indx1)
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
                  vec_c(indx2) = vechlp_c(indx1)
               END DO
            END DO
         END DO
      END DO
   end subroutine reorder_back_cmplx
END MODULE
