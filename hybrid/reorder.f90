MODULE m_reorder

CONTAINS
   SUBROUTINE reorder(nbasm, atoms, lcutm, maxlcutm, nindxm, imode, vec_r, vec_c)
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_atoms), INTENT(IN)   :: atoms

      ! - scalars -
      INTEGER, INTENT(IN)   ::  maxlcutm
      INTEGER, INTENT(IN)   ::  nbasm
      INTEGER, INTENT(IN)   ::  imode

      ! - arrays -
      INTEGER, INTENT(IN)   ::  lcutm(:)
      INTEGER, INTENT(IN)   ::  nindxm(0:maxlcutm, atoms%ntype)
      REAL, INTENT(INOUT), OPTIONAL::  vec_r(nbasm)
      COMPLEX, INTENT(INOUT), OPTIONAL::  vec_c(nbasm)
      ! - local scalars -
      INTEGER               ::  itype, ieq
      INTEGER               ::  indx1, indx2
      INTEGER               ::  l
      INTEGER               ::  n, m
      LOGICAL               :: l_real
      ! - local arrays -
      REAL                  ::  vechlp_r(nbasm)
      COMPLEX               ::  vechlp_c(nbasm)

      call timestart("reorder")
      l_real = PRESENT(vec_r)

      IF (imode /= 1 .and. imode /= 2) call judft_error('reorder: imode equals neither 1 nor 2')

      if (l_real) THEN
         vechlp_r = vec_r
      else
         vechlp_c = vec_c
      end if

      IF (imode == 1) THEN
         indx1 = 0
         indx2 = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               DO l = 0, lcutm(itype)
                  DO m = -l, l
                     DO n = 1, nindxm(l, itype) - 1
                        indx1 = indx1 + 1
                        indx2 = indx2 + 1
                        if (l_real) THEN
                           vec_r(indx1) = vechlp_r(indx2)
                        else
                           vec_c(indx1) = vechlp_c(indx2)
                        endif
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
                     if (l_real) THEN
                        vec_r(indx1) = vechlp_r(indx2)
                     else
                        vec_c(indx1) = vechlp_c(indx2)
                     endif

                  END DO
               END DO
            END DO
         END DO
      ELSE IF (imode == 2) THEN
         indx1 = 0
         indx2 = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               DO l = 0, lcutm(itype)
                  DO m = -l, l
                     DO n = 1, nindxm(l, itype) - 1
                        indx1 = indx1 + 1
                        indx2 = indx2 + 1
                        if (l_real) THEN
                           vec_r(indx2) = vechlp_r(indx1)
                        else
                           vec_c(indx2) = vechlp_c(indx1)
                        endif
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
                     if (l_real) THEN
                        vec_r(indx2) = vechlp_r(indx1)
                     else
                        vec_c(indx2) = vechlp_c(indx1)
                     endif
                  END DO
               END DO
            END DO
         END DO
      END IF
      !IR must not be rearranged
      call timestop("reorder")
   END SUBROUTINE reorder

END MODULE
