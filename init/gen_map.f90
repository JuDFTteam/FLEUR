MODULE m_gen_map
   IMPLICIT NONE

CONTAINS

!
! build up field map(iatom,isym), which contains the number of the atom, on
! which the atom iatom is mapped via the symmetry operation isym
! tvec is the translation, which maps R R_a + tau back in the unit cell
!
   SUBROUTINE gen_map(atoms, sym, oneD, hybinp)
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_oneD), INTENT(IN)  :: oneD
      TYPE(t_hybinp), INTENT(INOUT)::hybinp
      ! private scalars
      INTEGER                           :: iatom, first_eq_atom, itype, ieq, isym, iisym, ieq1
      INTEGER                           :: ratom, ok
      ! private arrays
      REAL                              :: rtaual(3)

      ALLOCATE (hybinp%map(atoms%nat, sym%nsym), stat=ok)
      IF (ok /= 0) call judft_error('gen_map: error during allocation of map')

      ALLOCATE (hybinp%tvec(3, atoms%nat, sym%nsym), stat=ok)
      IF (ok /= 0) call judft_error('gen_map: error during allocation of tvec')

      iatom = 0
      first_eq_atom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO isym = 1, sym%nsym

               IF (isym <= sym%nop) THEN
                  iisym = isym
               ELSE
                  iisym = isym - sym%nop
               END IF

               rtaual(:) = matmul(sym%mrot(:, :, iisym), atoms%taual(:, iatom)) + sym%tau(:, iisym)

               ratom = 0
               DO ieq1 = 1, atoms%neq(itype)
                  IF (all(abs(modulo(rtaual - atoms%taual(:, first_eq_atom + ieq1) + 1e-12, 1.0)) < 1e-10)) THEN
                     ratom = first_eq_atom + ieq1
                     hybinp%map(iatom, isym) = ratom
                     hybinp%tvec(:, iatom, isym) = nint(rtaual - atoms%taual(:, ratom))
                  END IF
               END DO
               IF (ratom == 0) call judft_error('eigen_hf: ratom not found')

            END DO
         END DO
         first_eq_atom = first_eq_atom + atoms%neq(itype)
      END DO

   END SUBROUTINE gen_map

END MODULE m_gen_map
