!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_jdos
   use m_judft
   use m_types_eigdos
   implicit none
   PRIVATE
   public t_jdos
   TYPE, extends(t_eigdos):: t_jDOS

      REAL, ALLOCATABLE    :: comp(:, :, :, :, :)  !decomposition in percent
      REAL, ALLOCATABLE    :: qmtp(:, :, :)      !How much of the state is in the muffin-tin sphere
      REAL, ALLOCATABLE    :: occ(:, :, :)       !Occupation of the j-states
      INTEGER, ALLOCATABLE  :: n_dos_to_na(:)

   CONTAINS
      PROCEDURE, PASS :: init => jDOS_init
      PROCEDURE      :: get_weight_eig
      PROCEDURE      :: get_num_weights
      PROCEDURE      :: get_weight_name
      PROCEDURE      :: get_spins
      PROCEDURE      :: calc_jDOS
      procedure      :: postprocessing
   END TYPE t_jDOS
CONTAINS
   subroutine postprocessing(this, noco,nococonv, banddos, alldos, ef)
      use m_types_atoms
      use m_types_noco
      use m_types_nococonv
      use m_types_banddos
      class(t_jDOS), intent(inout):: this
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_banddos), INTENT(IN)    :: banddos
      class(t_eigdos_list), intent(in), optional :: alldos(:)
      real, intent(in), optional :: ef
      return !currently no postprocessing needed for jdos
   end subroutine postprocessing    

   SUBROUTINE calc_jDOS(jDOS, ikpt, noccbd, ev_list, we, atoms, banddos, input, radfun, abc_u, abc_d)
      use m_types_atoms
      use m_types_banddos
      use m_types_input
      use m_types_radfun
      use m_types_abc
      use m_clebsch
      CLASS(t_jDOS), INTENT(INOUT)  :: jDOS
      TYPE(t_atoms), INTENT(IN)     :: atoms
      TYPE(t_banddos), INTENT(IN)     :: banddos
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_radfun), INTENT(IN)     :: radfun
      TYPE(t_abc), INTENT(IN)     :: abc_u, abc_d
      INTEGER, INTENT(IN)     :: ikpt
      INTEGER, INTENT(IN)     :: noccbd
      INTEGER, INTENT(IN)     :: ev_list(:)
      REAL, INTENT(IN)     :: we(:)

      INTEGER, PARAMETER :: lmax = 3 !Maximum l considered in j decomposition

      INTEGER :: n_dos, jcof, icof
      INTEGER :: jjj, iType, iBand, nn, iAtom, l, jj, j_ind, lmup, lmdown, spin, ilo, ilop
      REAL    :: j,mj, mup, mdown
      REAL    :: facup, facdown, summed, cf
      COMPLEX :: aup, bup, cup, adown, bdown, cdown, cupp, cdownp
      REAL    :: c(0:lmax*2)
      if (.not. jDOS%l_initialized) return
      DO iAtom = 1, atoms%nat
         iType = atoms%itype(iAtom)
         if (.not. banddos%dos_atom(iAtom)) cycle
         !find index for dos
         DO n_dos = 1, size(banddos%dos_atomlist)
            if (banddos%dos_atomlist(n_dos) == iAtom) exit
         END DO
         DO iBand = 1, noccbd
            j_ind = 0
            c = 0.0
            DO l = 0, lmax
            IF (l == 0) THEN
               !s-states (are not split up by SOC)
               DO jjj = 1, radfun%n_r(l)
                  DO jj = 1, radfun%n_r(l)
                   c(0) = c(0) + abc_u%cof(iband, 0, jjj, iatom)*conjg(abc_u%cof(iband, 0, jj, iatom))*radfun%integral(jjj, jj, 0, 1, 1)
                   c(0) = c(0) + abc_d%cof(iband, 0, jjj, iatom)*conjg(abc_d%cof(iband, 0, jj, iatom))*radfun%integral(jjj, jj, 0, 2, 2)

                  end do
               end do
            ELSE
               DO jj = 1, 2
                  j_ind = j_ind + 1
                  ! j = l +- 1/2
                  j = l + (jj - 1.5)
                  mj = -j
                  DO WHILE (mj <= j)
                     !mj = -l-+1/2, .... , l+-1/2

                     mup = mj - 0.5
                     mdown = mj + 0.5
                     DO icof = 1, radfun%n_r(l)
                        DO jcof = 1, radfun%n_r(l)

                           IF (input%jspins .EQ. 1) THEN
                              mdown = mdown*(-1)
                              spin = 1
                           ELSE
                              spin = 2
                           END IF

                           IF (ABS(mup) <= l) THEN
                              lmup = l*(l + 1) + INT(mup)
                              facup = clebsch(REAL(l), 0.5, mup, 0.5, j, mj)
                              aup = facup*abc_u%cof(iBand, lmup, icof, iAtom)
                              bup = facup*abc_u%cof(iBand, lmup, jcof, iAtom)
                           ELSE
                              aup = 0.0
                              bup = 0.0
                           END IF

                           IF (ABS(mdown) <= l) THEN
                              lmdown = l*(l + 1) + INT(mdown)
                              facdown = clebsch(REAL(l), 0.5, mdown, -0.5, j, mj)
                              adown = -1*facdown*abc_d%cof(iBand, lmdown, icof, iAtom)
                              bdown = -1*facdown*abc_d%cof(iBand, lmdown, jcof, iAtom)
                           ELSE
                              adown = 0.0
                              bdown = 0.0
                           END IF

                           !c := norm of facup |up> + facdown |down>
                           !We have to write it out explicitely because
                           !of the offdiagonal scalar products that appear
                           c(j_ind) = c(j_ind) &
                                      + aup*CONJG(bup)*radfun%integral(icof, jcof, l, 1, 1) &
                                      + adown*CONJG(bdown)*radfun%integral(icof, jcof, l, 2, 2) &
                                      + aup*CONJG(bdown)*radfun%integral(icof, jcof, l, 1, 2) &
                                      + adown*CONJG(bup)*radfun%integral(icof, jcof, l, 2, 1)
                        end do
                     end do
                     
                     mj = mj + 1
                  END DO
               END DO
            END IF
            END DO
            summed = SUM(c(0:2*lmax))
            cf = 100.0/summed
            j_ind = 0
            DO l = 0, 3
            DO jj = 1, 2
               IF (l /= 0) j_ind = j_ind + 1
               jDOS%comp(ev_list(iBand), l, jj, n_dos, ikpt) = c(j_ind)*cf
               jDOS%qmtp(ev_list(iBand), n_dos, ikpt) = 100.0*summed
               jDOS%occ(l, jj, iAtom) = jDOS%occ(l, jj, n_dos) + we(iBand)*c(j_ind)
            END DO
            END DO
         END DO
      END DO

   END SUBROUTINE calc_jDOS

   pure integer function get_spins(this)
      CLASS(t_jdos), INTENT(IN)::this
      get_spins = 1
   END function

   function get_weight_eig(this, id)
      class(t_jdos), intent(in):: this
      INTEGER, intent(in)      :: id
      real, allocatable:: get_weight_eig(:, :, :)

      integer :: i, l, jj, na

      ALLOCATE (get_weight_eig(size(this%comp, 1), size(this%comp, 5), 1))

      i = 0
      DO na = 1, size(this%comp, 4)
         DO l = 0, 3
            DO jj = 1, MERGE(1, 2, l == 0)
               i = i + 1
               if (i == id) get_weight_eig(:, :, 1) = this%comp(:, l, jj, na, :)*this%qmtp(:, na, :)/10000.
               if (i > id) RETURN
            END DO
         END DO
      END DO
   end function

   integer function get_num_weights(this)
      class(t_jdos), intent(in):: this
      get_num_weights = 7*size(this%comp, 4)
   end function

   character(len=20) function get_weight_name(this, id)
      class(t_jdos), intent(in):: this
      INTEGER, intent(in)         :: id
      integer :: i, l, jj, na
      character :: spdfg(0:4) = ["s", "p", "d", "f", "g"]
      character(len=3) :: jname

      i = 0
      DO na = 1, size(this%comp, 4)
         DO l = 0, 3
            DO jj = -1, MERGE(-1, 1, l == 0), 2
               i = i + 1
               WRITE (jname, '(i1,a,i1)') INT(2*l + jj), '-', 2
               if (i == id) THEN
                  IF (l .EQ. 0) write (get_weight_name, "(a,i0,a)") "jDOS:", this%n_dos_to_na(na), spdfg(l)
                  IF (l .NE. 0) write (get_weight_name, "(a,i0,a,a)") "jDOS:", this%n_dos_to_na(na), spdfg(l), jname
               end if
               if (i > id) RETURN
            END DO
         END DO
      END DO

   end function

   SUBROUTINE jDOS_init(thisjDOS, input, banddos, atoms, kpts, eig)

      USE m_types_setup
      USE m_types_kpts

      IMPLICIT NONE

      CLASS(t_jDOS), INTENT(INOUT) :: thisjDOS
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_banddos), INTENT(IN)    :: banddos

      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_kpts), INTENT(IN)    :: kpts
      REAL, INTENT(IN)                      :: eig(:, :, :)
      thisjDOS%l_initialized = .TRUE.
      thisjDOS%n_dos_to_na = banddos%dos_atomlist
      IF (banddos%l_jdos .AND. banddos%dos) THEN
         ALLOCATE (thisjDOS%comp(input%neig, 0:3, 2, size(banddos%dos_atomlist), kpts%nkpt), source=0.0)
         ALLOCATE (thisjDOS%qmtp(input%neig, size(banddos%dos_atomlist), kpts%nkpt), source=0.0)
         ALLOCATE (thisjDOS%occ(0:3, 2, size(banddos%dos_atomlist)), source=0.0)
         thisjDOS%eig = eig
      ELSE
         ALLOCATE (thisjDOS%dos(0, 0, 0))
         ALLOCATE (thisjDOS%comp(1, 1, 1, 0, 1), source=0.0)
         ALLOCATE (thisjDOS%qmtp(1, 0, 1), source=0.0)
         ALLOCATE (thisjDOS%occ(1, 1, 0), source=0.0)
      END IF

      thisjDOS%name_of_dos = "jDOS"

   END SUBROUTINE jDOS_init
end module m_types_jDOS
