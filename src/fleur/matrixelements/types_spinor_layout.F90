!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_spinor_layout
   !> How the eigenvectors of one k-point are laid out, and how many radial sets
   !> exist to contract them with.
   !>
   !> The number of records in the eig file does not say what the records are: two
   !> records are the two halves of one spinor after a second variation, but two
   !> independent spin channels in a collinear calculation, and stacking the latter
   !> is meaningless. One record is a scalar state without spin structure, but a
   !> whole spinor of 2N rows when the Hamiltonian was set up non-collinearly. This
   !> type carries the distinction next to the count, so a consumer asks instead of
   !> re-deriving it from the input flags.
   !>
   !> radial_slot is a function of the radial array rather than a method, because the
   !> number of radial sets is a property of that array and every caller already holds
   !> it. Nothing else about the layout is needed to answer it.

   USE m_judft
   USE m_types_input
   USE m_types_noco
   USE m_types_lapw
   USE m_types_atoms
   USE m_types_radfun

   IMPLICIT NONE
   PRIVATE

   !> What the records mean. Distinct from how many there are.
   INTEGER, PARAMETER, PUBLIC :: LAYOUT_SCALAR   = 1  !> one state, no spin index
   INTEGER, PARAMETER, PUBLIC :: LAYOUT_CHANNELS = 2  !> two independent spin channels
   INTEGER, PARAMETER, PUBLIC :: LAYOUT_SPINOR   = 3  !> two components of one state

   TYPE t_spinor_layout
      INTEGER :: layout    = -1       !> LAYOUT_SCALAR, LAYOUT_CHANNELS or LAYOUT_SPINOR
      INTEGER :: nrec      = -1       !> records to read from the eig file: 1 or 2
      LOGICAL :: l_stacked = .FALSE.  !> a spinor that already arrives as one 2N matrix
      INTEGER :: row_dn    = 0        !> first spin-down row of a stacked spinor, else 0
      INTEGER :: n_radial  = -1       !> radial sets available to contract with
   CONTAINS
      PROCEDURE :: init       => spinor_layout_init
      PROCEDURE :: nbasfcn    => spinor_layout_nbasfcn
      PROCEDURE :: row_offset => spinor_layout_row_offset
   END TYPE t_spinor_layout

   PUBLIC :: t_spinor_layout, radial_slot, n_radial_of

CONTAINS

   SUBROUTINE spinor_layout_init(this, input, noco, lapw, atoms, radfun, l_both_spinors)
      !> l_both_spinors distinguishes the two eig files a spin-orbit run produces:
      !> before the second variation it holds input%jspins independent channels,
      !> after it two spinor records. That is a property of the stage rather than
      !> of any input flag, so the caller has to state it.
      CLASS(t_spinor_layout), INTENT(OUT) :: this
      TYPE(t_input),          INTENT(IN)  :: input
      TYPE(t_noco),           INTENT(IN)  :: noco
      TYPE(t_lapw),           INTENT(IN)  :: lapw
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      !> The radial functions the contraction will index. Their fourth extent is the
      !> only honest source of truth about how many sets exist: input%jspins says how
      !> many were asked for, which is not the same statement.
      TYPE(t_radfun), OPTIONAL, INTENT(IN) :: radfun(:)
      LOGICAL,        OPTIONAL, INTENT(IN) :: l_both_spinors

      LOGICAL :: l_spinor_records

      l_spinor_records = .FALSE.
      IF (PRESENT(l_both_spinors)) l_spinor_records = l_both_spinors

      IF (noco%l_noco) THEN
         !> The Hamiltonian was 2N x 2N, so one record already holds the whole spinor.
         this%layout    = LAYOUT_SPINOR
         this%nrec      = 1
         this%l_stacked = .TRUE.
      ELSE IF (l_spinor_records) THEN
         !> Second variation done: two records, and they are halves of one state.
         this%layout    = LAYOUT_SPINOR
         this%nrec      = 2
         this%l_stacked = .FALSE.
      ELSE IF (input%jspins == 2) THEN
         !> Two potentials, two eigenproblems. The records are unrelated states.
         this%layout    = LAYOUT_CHANNELS
         this%nrec      = 2
         this%l_stacked = .FALSE.
      ELSE
         this%layout    = LAYOUT_SCALAR
         this%nrec      = 1
         this%l_stacked = .FALSE.
      END IF

      !> Only a spinor has a spin-down block to point at.
      this%row_dn = 0
      IF (this%layout == LAYOUT_SPINOR) this%row_dn = lapw%nv(1) + atoms%nlotot

      IF (PRESENT(radfun)) THEN
         this%n_radial = n_radial_of(radfun)
      ELSE
         this%n_radial = input%jspins
      END IF

      IF (l_spinor_records .AND. (noco%l_noco .OR. .NOT.noco%l_soc)) &
         CALL judft_error("t_spinor_layout: spinor records only exist for l_soc=T, l_noco=F", &
                          calledby="spinor_layout_init")
      IF (this%row_dn == 0 .AND. this%layout == LAYOUT_SPINOR) &
         CALL judft_error("t_spinor_layout: a spinor needs a spin-down row", &
                          calledby="spinor_layout_init")
      IF (this%n_radial < 1 .OR. this%n_radial > 2) &
         CALL judft_error("t_spinor_layout: there are neither one nor two radial sets", &
                          calledby="spinor_layout_init")
      IF (this%layout == LAYOUT_CHANNELS .AND. this%n_radial /= 2) &
         CALL judft_error("t_spinor_layout: two spin channels without two radial sets", &
                          calledby="spinor_layout_init")
   END SUBROUTINE spinor_layout_init


   PURE INTEGER FUNCTION n_radial_of(radfun)
      !> How many radial sets exist, read from the array that will be indexed. This is
      !> not the same statement as input%jspins, which says how many were asked for, and
      !> the point of reading it here is that the two can disagree.
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      n_radial_of = SIZE(radfun(1)%integral, 4)
   END FUNCTION n_radial_of


   PURE INTEGER FUNCTION radial_slot(radfun, isp)
      !> The radial set to contract spin component isp with. Clamped, because a spinor
      !> has two components while a single potential generates one radial set: asking
      !> for slot 2 there reads past the array.
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      INTEGER,        INTENT(IN) :: isp
      radial_slot = MIN(isp, n_radial_of(radfun))
   END FUNCTION radial_slot


   PURE INTEGER FUNCTION spinor_layout_row_offset(this, isp)
      !> Row at which spin component isp begins inside one record. Zero unless the
      !> spinor arrives stacked, where the second component starts at row_dn.
      CLASS(t_spinor_layout), INTENT(IN) :: this
      INTEGER,                INTENT(IN) :: isp
      spinor_layout_row_offset = 0
      IF (this%l_stacked .AND. isp == 2) spinor_layout_row_offset = this%row_dn
   END FUNCTION spinor_layout_row_offset


   PURE INTEGER FUNCTION spinor_layout_nbasfcn(this, lapw, atoms)
      !> Rows of one eigenvector record: both spin blocks when the spinor arrives
      !> stacked, one otherwise.
      CLASS(t_spinor_layout), INTENT(IN) :: this
      TYPE(t_lapw),           INTENT(IN) :: lapw
      TYPE(t_atoms),          INTENT(IN) :: atoms
      IF (this%l_stacked) THEN
         spinor_layout_nbasfcn = lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot
      ELSE
         spinor_layout_nbasfcn = lapw%nv(1) + atoms%nlotot
      END IF
   END FUNCTION spinor_layout_nbasfcn

END MODULE m_types_spinor_layout
