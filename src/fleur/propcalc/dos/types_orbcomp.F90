!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_orbcomp
   use m_judft
   use m_types_eigdos
   implicit none
   PRIVATE
   integer, allocatable:: combine(:, :)
   real, allocatable   :: factor(:, :)
   public t_orbcomp
   TYPE, extends(t_eigdos):: t_orbcomp

      REAL, ALLOCATABLE    :: comp(:, :, :, :, :)
      REAL, ALLOCATABLE    :: qmtp(:, :, :, :)
      INTEGER, ALLOCATABLE  :: n_dos_to_na(:)
   CONTAINS
      PROCEDURE, PASS :: init => orbcomp_init
      PROCEDURE      :: get_num_weights
      PROCEDURE      :: get_weight_eig
      PROCEDURE      :: get_weight_name
      PROCEDURE      :: sym_weights
      PROCEDURE      :: calc_orb_comp
      PROCEDURE      :: postprocessing
   END TYPE t_orbcomp
CONTAINS

 subroutine postprocessing(this, noco,nococonv, banddos)
      use m_types_atoms
      use m_types_noco
      use m_types_nococonv
      use m_types_banddos
      class(t_orbcomp), intent(inout):: this
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_banddos), INTENT(IN)    :: banddos
      return !currently no postprocessing needed for orbcomp
   end subroutine postprocessing 
   subroutine set_coefficients()
    !! Here the coefficients and factor for the linear combinations are defined that
    !! generate the correct orbital resolved DOS
    !! Up to four different coefficients might be combined and there exist 16 different lm up to the f-states.
      real, parameter::h = sqrt(0.5), g = sqrt(0.0625)
      if (allocated(combine)) return !coefficients present
      allocate (combine(4, 23), factor(4, 23))
      combine = -1; factor = 1
!s-state
      combine(1, 1) = 0
!p-states
      combine(1, 2) = 1; combine(2, 2) = 3
      factor(1, 2) = h; factor(2, 2) = -h
      combine(1, 3) = 1; combine(2, 3) = 3
      factor(1, 3) = h; factor(2, 3) = h
      combine(1, 4) = 2
!d-states
      combine(1, 5) = 4; combine(2, 5) = 8
      factor(1, 5) = h; factor(2, 5) = -h
      combine(1, 6) = 5; combine(2, 6) = 7
      factor(1, 6) = h; factor(2, 6) = h
      combine(1, 7) = 5; combine(2, 7) = 7
      factor(1, 7) = h; factor(2, 7) = -h
      combine(1, 8) = 4; combine(2, 8) = 8
      factor(1, 8) = h; factor(2, 8) = h
      combine(1, 9) = 6
!f-states: a cubic set (cub)
      combine(1, 10) = 9; combine(2, 10) = 15; combine(3, 10) = 11; combine(4, 10) = 13
      factor(1, 10) = sqrt(5.)*g; factor(2, 10) = -sqrt(5.)*g; factor(3, 10) = -sqrt(3.0)*g; factor(4, 10) = sqrt(3.0)*g
      combine(1, 11) = 9; combine(2, 11) = 15; combine(3, 11) = 11; combine(4, 11) = 13
      factor(1, 11) = sqrt(5.)*g; factor(2, 11) = sqrt(5.)*g; factor(3, 11) = sqrt(3.0)*g; factor(4, 11) = sqrt(3.0)*g
      combine(1, 12) = 12
      combine(1, 13) = 9; combine(2, 13) = 15; combine(3, 13) = 11; combine(4, 13) = 13
      factor(1, 13) = sqrt(3.)*g; factor(2, 13) = sqrt(3.)*g; factor(3, 13) = -sqrt(5.0)*g; factor(4, 13) = -sqrt(5.0)*g
      combine(1, 14) = 10; combine(2, 14) = 14
      factor(1, 14) = h; factor(2, 14) = h
      combine(1, 15) = 9; combine(2, 15) = 15; combine(3, 15) = 11; combine(4, 15) = 13
      factor(1, 15) = sqrt(3.)*g; factor(2, 15) = -sqrt(3.)*g; factor(3, 15) = sqrt(5.0)*g; factor(4, 15) = -sqrt(5.0)*g
      combine(1, 16) = 10; combine(2, 16) = 14
      factor(1, 16) = h; factor(2, 16) = -h
!f-states:        a low symmetry set (lss) (not used by default)
      combine(1, 17) = 11; combine(2, 17) = 13
      factor(1, 17) = h; factor(2, 17) = -h
      combine(1, 18) = 11; combine(2, 18) = 13
      factor(1, 18) = h; factor(2, 18) = h
      combine(1, 19) = 12
      combine(1, 20) = 10; combine(2, 20) = 14
      factor(1, 20) = h; factor(2, 20) = -h
      combine(1, 21) = 10; combine(2, 21) = 14
      factor(1, 21) = h; factor(2, 21) = h
      combine(1, 22) = 9; combine(2, 22) = 15
      factor(1, 22) = h; factor(2, 22) = -h
      combine(1, 23) = 9; combine(2, 23) = 15
      factor(1, 23) = h; factor(2, 23) = h
   end subroutine

   SUBROUTINE calc_orb_comp(orbcomp, atoms, banddos, radfun, abc_in, abc1_in, ev_list, itype, ikpt, jsp, jsp1)
    !!     Calculates an orbital composition of eigen states
    !! Based on code from    Yury  Koroteev  2003-12-24

      USE m_types_atoms
      USE m_types_banddos
      USE m_types_abc
      USE m_types_radfun

      IMPLICIT NONE
      CLASS(t_orbcomp), INTENT(INOUT)  :: orbcomp
      TYPE(t_atoms), INTENT(IN)        :: atoms
      TYPE(t_banddos), INTENT(IN)      :: banddos
      TYPE(t_radfun), INTENT(IN)      :: radfun
      Type(t_abc), intent(in), target    :: abc_in, abc1_in
      integer, intent(in)               :: itype, ikpt, jsp, jsp1
      INTEGER, intent(in)                  :: ev_list(:)

      COMPLEX:: comp(23)
      type(t_abc), TARGET :: abc_rot, abc1_rot
      type(t_abc), POINTER::abc, abc1

      INTEGER:: n, i, n_orb, l, j, jj, n_dos, mt, natom, jspin
      complex:: sc, lin_comb, lin_comb1

      call set_coefficients()
      jspin = merge(jsp, 3, jsp == jsp1)

      do natom = 1, atoms%neq(itype)
         mt = natom + atoms%firstatom(itype) - 1
         !find index for dos
         DO n_dos = 1, size(banddos%dos_atomlist)
            if (banddos%dos_atomlist(n_dos) == mt) exit
         END DO
         if (n_dos > size(banddos%dos_atomlist)) cycle ! no n_dos for this atom found
         IF (ANY((/banddos%alpha(mt), banddos%beta(mt), banddos%gamma(mt)/) .NE. 0.0)) THEN !check if atom should be rotated....
            abc_rot=abc%rotate(banddos%alpha(mt), banddos%beta(mt), banddos%gamma(mt),3)
            abc1_rot=abc1%rotate(banddos%alpha(mt), banddos%beta(mt), banddos%gamma(mt),3)
            abc => abc_rot
            abc1 => abc1_rot
            !CALL abcrot2(ityp, mt, atoms, banddos, eigVecCoeffs, jspin, acof, bcof, ccof) ! rotate ab-coeffs
         ELSE
            abc => abc_in
            abc1 => abc1_in
         end if
         DO i = 1, size(abc%cof, 1) !loop over all energies
            comp = CMPLX(0., 0.)
            DO n_orb = 1, size(combine, 2) !loop over all orbital decompositions
               l = floor(sqrt(1.*combine(1, n_orb)))
               DO j = 1, size(abc%cof, 3)
                  DO jj = 1, size(abc%cof, 3)
                     lin_comb=0.0;lin_comb1=0.0
                     !now create linear combination of abc coefficients
                     DO n = 1, size(combine, 1)
                        if (combine(n, n_orb) == -1) cycle !no more coeffients to combine
                        lin_comb = lin_comb + factor(n, n_orb)*abc%cof(i, combine(n, n_orb), j, natom)
                        lin_comb1 = lin_comb1 + factor(n, n_orb)*abc1%cof(i, combine(n, n_orb), jj, natom)
                     end do
                     comp(n_orb) = comp(n_orb) + lin_comb*CONJG(lin_comb1)*radfun%integral(j, jj, l, jsp, jsp1)
                  END DO
               END DO
            END DO

            !    calculate an orbital cnomposition in percets
            !    Only the first 16 decompositions are used to create the sum
            sc = sum(comp(1:16))
            orbcomp%qmtp(ev_list(i), n_dos, ikpt, jspin) = 100*sc
            if (jspin == 3) orbcomp%qmtp(ev_list(i), n_dos, ikpt, 4) = aimag(100*sc)
            if (abs(sc) > 1E-18) then
               orbcomp%comp(ev_list(i), :, n_dos, ikpt, jspin) = comp(:)*100./sc
               if (jspin == 3) orbcomp%comp(ev_list(i), :, n_dos, ikpt, 4) = aimag(comp(:)*100./sc)
            end if
         END DO !loop over energies
      end do!atom loop

   END SUBROUTINE calc_orb_comp

   subroutine sym_weights(this)
      class(t_orbcomp), intent(inout):: this
      integer:: i, j
      return !This is done later in get_weights
   end subroutine

   integer function get_num_weights(this)
      class(t_orbcomp), intent(in):: this
      get_num_weights = 23*size(this%comp, 3)
   END function

   character(len=20) function get_weight_name(this, id)
      class(t_orbcomp), intent(in):: this
      INTEGER, intent(in)         :: id

      INTEGER :: ind, na, nc
      ind = 0
      DO na = 1, size(this%comp, 3)
         DO nc = 1, 23
            ind = ind + 1
            if (ind == id) THEN
               write (get_weight_name, "(a,i0,a,i0)") "ORB:", this%n_dos_to_na(na), ",ind:", nc
               RETURN
            ELSE IF (ind > id) then
               CALL judft_bug("Types_orbcomp: data not found")
            END IF
         END DO
      END DO
   end function

   function get_weight_eig(this, id)
      class(t_orbcomp), intent(in):: this
      INTEGER, intent(in)      :: id
      real, allocatable:: get_weight_eig(:, :, :)

      integer :: i, ind, na

      ind = 0
      DO na = 1, size(this%comp, 3)
         DO i = 1, 23
            ind = ind + 1
            if (ind == id) THEN
               get_weight_eig = this%comp(:, i, na, :, :)*this%qmtp(:, na, :, :)/10000.
               call this%sym_weights_eigdos(get_weight_eig)
               return
            END IF
         END DO
      END DO
   end function

   SUBROUTINE orbcomp_init(thisOrbcomp, input, banddos, atoms, kpts, eig)

      USE m_types_setup
      USE m_types_kpts

      IMPLICIT NONE

      CLASS(t_orbcomp), INTENT(INOUT) :: thisOrbcomp
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_banddos), INTENT(IN)    :: banddos

      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_kpts), INTENT(IN)    :: kpts
      REAL, INTENT(IN)                      :: eig(:, :, :)
      thisOrbcomp%n_dos_to_na = banddos%dos_atomlist
      IF ((banddos%l_orb) .AND. banddos%dos) THEN
         ALLOCATE (thisOrbcomp%comp(input%neig, 23, size(banddos%dos_atomlist), kpts%nkpt, input%jspins))
         ALLOCATE (thisOrbcomp%qmtp(input%neig, size(banddos%dos_atomlist), kpts%nkpt, input%jspins))
         thisOrbcomp%eig = eig
      ELSE
         ALLOCATE (thisOrbcomp%dos(0, 0, 0))
         ALLOCATE (thisOrbcomp%comp(1, 1, 0, 1, input%jspins))
         ALLOCATE (thisOrbcomp%qmtp(1, 0, 1, input%jspins))
      END IF

      thisOrbcomp%comp = 0.0
      thisOrbcomp%qmtp = 0.0
      thisOrbcomp%name_of_dos = "Orbcomp"
   END SUBROUTINE orbcomp_init

end module m_types_orbcomp
