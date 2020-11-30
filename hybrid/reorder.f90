MODULE m_reorder
   interface reorder_forw 
      module procedure reorder_forw_real, reorder_forw_cmplx
   end interface reorder_forw

   interface reorder_back
      module procedure reorder_back_real, reorder_back_cmplx
   end interface reorder_back

CONTAINS
   subroutine reorder_forw_real(nbasm, atoms, lcutm, nindxm, vec_r)
      USE m_types
      USE m_juDFT
      IMPLICIT NONE

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
                     vec_r(indx1) = vechlp_r(indx2)
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
                  vec_r(indx1) = vechlp_r(indx2)
               END DO
            END DO
         END DO
      END DO
   end subroutine reorder_forw_real

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

   subroutine reorder_forw_cmplx(nbasm, atoms, lcutm, nindxm, vec_c)
      USE m_types
      USE m_juDFT
      use m_constants, only: cmplx_0
      IMPLICIT NONE

      INTEGER, INTENT(IN)       :: nbasm, lcutm(:), nindxm(0:, :)
      TYPE(t_atoms), INTENT(IN) :: atoms
      complex, INTENT(INOUT)    :: vec_c(nbasm)
      
      INTEGER               :: itype, ieq, indx1, indx2, l, n, m, info
      complex, allocatable  ::  vechlp_c(:)

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
                     vec_c(indx1) = vechlp_c(indx2)
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
                  vec_c(indx1) = vechlp_c(indx2)
               END DO
            END DO
         END DO
      END DO
   end subroutine reorder_forw_cmplx

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
