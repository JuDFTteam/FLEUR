!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_dos
   USE m_juDFT
   USE m_types_eigdos
   IMPLICIT NONE
   PRIVATE
   PUBLIC:: t_dos
   TYPE, extends(t_eigdos):: t_dos
      INTEGER, ALLOCATABLE :: neq(:)
      INTEGER, ALLOCATABLE :: jsym(:, :, :)
      REAL, ALLOCATABLE :: qis(:, :, :)
      REAL, ALLOCATABLE :: qal(:, :, :, :, :)
      REAL, ALLOCATABLE :: qTot(:, :, :)
      CHARACTER(len=20), ALLOCATABLE:: weight_names(:)!This must be allocated in init of derived type

   CONTAINS
      PROCEDURE, PASS :: init => dos_init
      PROCEDURE      :: get_weight_eig
      PROCEDURE      :: get_num_spins
      PROCEDURE      :: get_num_weights
      PROCEDURE      :: get_weight_name
      PROCEDURE      :: sym_weights
      PROCEDURE      :: calc_mt_dos
      PROCEDURE      :: postprocessing
   END TYPE t_dos

CONTAINS

   subroutine postprocessing(this, noco,nococonv, banddos)
      use m_types_atoms
      use m_types_noco
      use m_types_nococonv
      use m_types_banddos
      class(t_dos), intent(inout):: this
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_banddos), INTENT(IN)    :: banddos


      integer:: n_dos,ikpt,i,l
      complex:: qal21
      if (.not.noco%l_noco) return
      if (banddos%global_frame) THEN !Only if global frame is requested do the transformation
         DO n_dos = 1, size(this%qal, 2)
            if (banddos%dos_typelist(n_dos) == 0) cycle !this is not a  MT DOS
            do ikpt = 1, size(this%qal, 4)
               DO i = 1, size(this%qal, 3)
                  DO l = 0, 3
                     qal21 = cmplx(this%qal(l, n_dos, i, ikpt, 3), this%qal(l, n_dos, i, ikpt, 4))
                     CALL nococonv%rotdenmat(nococonv%alph(banddos%dos_typelist(n_dos)), nococonv%beta(banddos%dos_typelist(n_dos)), &
                                   this%qal(l, n_dos, i, ikpt, 1), this%qal(l, n_dos, i, ikpt, 2), qal21,.true.)
                     this%qal(l, n_dos, i, ikpt, 3) = real(qal21)
                     this%qal(l, n_dos, i, ikpt, 4) = aimag(qal21)
                  END DO
               END DO
            END DO
         END DO
      endif   
   end subroutine

   subroutine calc_mt_dos(dos, abc, abc1, banddos, radfun, atoms, ev_list, itype, ikpt, jsp, jsp1)
    use m_types_atoms
    use m_types_banddos
    use m_types_radfun
    use m_types_abc

    class(t_dos), intent(inout):: dos
    TYPE(t_atoms), INTENT(IN)    :: atoms
    TYPE(t_abc), INTENT(IN)    :: abc, abc1
    TYPE(t_banddos), INTENT(IN)    :: banddos
    TYPE(t_radfun), INTENT(IN)    :: radfun
    INTEGER, INTENT(IN) :: ev_list(:)
    INTEGER, INTENT(IN)   :: itype, ikpt, jsp, jsp1

    complex:: suma
    INTEGER:: i, l, lm, j, n_dos, m, natom,jj
    Integer:: jspin

    jspin = merge(jsp, 3, jsp == jsp1)
    DO i = 1, size(abc%cof, 1)
       DO l = 0, atoms%lmax(itype)
          suma = CMPLX(0., 0.)
          DO m = -l, l
             lm = l*(l + 1) + m
             DO natom = 1, atoms%neq(itype)
                DO j = 1, radfun%n_r(l)
                   DO jj = 1, radfun%n_r(l)
                     suma = suma + abc%cof(i, lm, j, natom)*CONJG(abc1%cof(i, lm, jj, natom))*radfun%integral(j, jj, l, jsp, jsp1)
                   END DO
                END DO
             END DO
          END DO
      
          if (l < 4) THEN !Only if l<3 store local projects
             DO n_dos = 1, size(banddos%dos_typelist)
                if (itype == banddos%dos_typelist(n_dos)) THEN !Only store projected DOS if this atom is in the dos_typelist
                   dos%qal(l, n_dos, ev_list(i), ikpt, jspin) = suma/atoms%neq(itype)
                   if (jspin == 3) dos%qal(l, n_dos, ev_list(i), ikpt, 4) = aimag(suma)/atoms%neq(itype)
                end if
             end do
          end if
          dos%qTot(ev_list(i), ikpt, jspin) = dos%qTot(ev_list(i), ikpt, jsp) + suma
          if (jspin == 3) dos%qTot(ev_list(i), ikpt, 4) = dos%qTot(ev_list(i), ikpt, 4) + aimag(suma)
       END DO
    END DO
   end subroutine

   subroutine sym_weights(this)
      class(t_dos), intent(inout):: this

      integer:: i, j

      DO i = 0, size(this%qal, 1) - 1
         DO j = 1, size(this%qal, 2)
            call this%sym_weights_eigdos(this%qal(i, j, :, :, :))
         end do
      END DO
      call this%sym_weights_eigdos(this%qis(:, :, :))
      call this%sym_weights_eigdos(this%qtot(:, :, :))
   end subroutine

   integer function get_num_weights(this)
      class(t_dos), intent(in):: this
      get_num_weights = 0
      if (allocated(this%weight_names)) get_num_weights = size(this%weight_names)
   end function

   character(len=20) function get_weight_name(this, id)
      class(t_dos), intent(in):: this
      INTEGER, intent(in)         :: id
      if (.not. allocated(this%weight_names)) call judft_error("No weight names in t_eigdos")
      if (id > size(this%weight_names)) call judft_error("Not enough weight names in t_eigdos")
      get_weight_name = this%weight_names(id)
   end function

   integer function get_num_spins(this)
      class(t_dos), intent(in):: this
      get_num_spins = size(this%qis, 3)
   end function

   function get_weight_eig(this, id)
      class(t_dos), intent(in):: this
      INTEGER, intent(in)     :: id
      real, allocatable:: get_weight_eig(:, :, :)

      INTEGER :: ind, l, ntype, i
      allocate (get_weight_eig, mold=this%qis)

      if (id == 1) THEN
         get_weight_eig = this%qTot
         if (all(this%qis == 0.0)) then
            get_weight_eig = 1.0
         END IF
      END IF
      if (id == 2) THEN
         get_weight_eig = this%qis
         if (all(get_weight_eig == 0.0)) then
            !No INT dos calculated so far...
            get_weight_eig = 1.0
            DO ntype = 1, size(this%qal, 2)
               DO l = 0, 3
                  get_weight_eig = get_weight_eig - this%qal(l, ntype, :, :, :)*this%neq(ntype)
               END DO
            END DO
         end if
      end if
      if (id == 3) get_weight_eig = 1.*this%jsym
      ind = 3
      DO ntype = 1, size(this%qal, 2)
         DO l = 0, 3
            ind = ind + 1
            if (ind == id) get_weight_eig = this%qal(l, ntype, :, :, :)
         END DO
      END DO
   end function

   SUBROUTINE dos_init(thisDOS, input, atoms, kpts, banddos, l_noco, eig)
      USE m_types_input
      USE m_types_atoms
      USE m_types_banddos
      USE m_types_kpts
      USE m_types_noco
      IMPLICIT NONE
      CLASS(t_dos), INTENT(INOUT) :: thisDOS
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_atoms), INTENT(IN)    :: atoms
      LOGICAL, INTENT(IN)          :: l_noco
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_banddos), INTENT(IN)    :: banddos
      real, intent(in)                       :: eig(:, :, :)

      INTEGER :: ntype, l, i, ind,ispin
      character :: spdfg(0:4) = ["s", "p", "d", "f", "g"]
      thisDOS%name_of_dos = "Local"
      thisDOS%neq = atoms%neq(banddos%dos_typelist)
      thisDOS%eig = eig
      ispin= merge(4,input%jspins,l_noco)
      ALLOCATE (thisDOS%jsym(input%neig, kpts%nkpt, ispin))
      ALLOCATE (thisDOS%qis(input%neig, kpts%nkpt, ispin))
      ALLOCATE (thisDOS%qal(0:3, size(banddos%dos_typelist), input%neig, kpts%nkpt,ispin))
      ALLOCATE (thisDOS%qTot(input%neig, kpts%nkpt, ispin))

      thisDOS%jsym = 0
      thisDOS%qis = 0.0
      thisDOS%qal = 0.0
      thisDOS%qTot = 0.0

      allocate (thisDOS%weight_names(3 + 4*size(banddos%dos_typelist)))
      thisDOS%weight_names(1) = "Total"
      thisDOS%weight_names(2) = "INT"
      thisDOS%weight_names(3) = "Sym"
      ind = 3
      DO ntype = 1, size(banddos%dos_typelist)
         DO l = 0, 3
            ind = ind + 1
            write (thisDOS%weight_names(ind), "(a,i0,a)") "MT:", banddos%dos_typelist(ntype), spdfg(l)
         END DO
      END DO

   END SUBROUTINE dos_init

END MODULE m_types_dos
