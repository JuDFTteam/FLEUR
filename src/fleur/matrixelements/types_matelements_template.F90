!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  TEMPLATE for a new matrix-element provider. Copy this file, rename it, and
!>  replace the angular factor; the rest is the contract and can stay as it is.
!>
!>  KEEP THIS FILE IN THE BUILD, although nothing calls it. A skeleton that is not
!>  compiled stops matching the interface it claims to demonstrate, and nothing says
!>  so -- which is how the schema came to accept two operator names the tables do not
!>  have. Compiled, a change to t_matelements breaks the build on the commit that
!>  makes it. It costs one object file, and its CMakeLists entry is a commit of its
!>  own so that dropping the pattern is a single revert.
!>
!>  As it stands it computes the MUFFIN-TIN OVERLAP, i.e. the identity restricted to
!>  the spheres:
!>
!>      S_mn = < psi_m | psi_n >_MT
!>
!>  which is a real quantity and one you can check by eye: its diagonal is the
!>  muffin-tin charge of the state, so it is real, positive and below one. Start from
!>  a kernel whose answer you already know, then put your own angular factor in.
!>
!>  The shape demonstrated here -- spin-diagonal, one site, three Cartesian
!>  components -- is the most common one. The live example of the same shape with a
!>  real operator in it is m_types_matelements_orbital; read that one next.
!>
!>  See also: the tutorial in August/refactor/tutorial_operadores.pdf, and the
!>  checklist in this directory's README.md for the other six files to touch.
MODULE m_types_matelements_template
   USE m_types_matelements
   USE m_types_mat
   USE m_types_abc
   USE m_types_radfun
   USE m_types_spinor_layout, ONLY: radial_slot
   USE m_types_usdus
   USE m_types_atoms
   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_matelements) :: t_matelements_template
      TYPE(t_atoms), POINTER :: atoms => NULL()
      INTEGER :: ntyp    = 0   !> atom type
      INTEGER :: iat     = 0   !> which of the equivalent atoms of that type
      INTEGER :: channel = 0   !> which spin channel this instance serves; 0 = a spinor
   CONTAINS
      PROCEDURE :: init                 => melem_template_init
      PROCEDURE :: calc_matrix_elements => melem_template_calc
   END TYPE t_matelements_template

   PUBLIC :: t_matelements_template

CONTAINS

   !> Bind the object to one atom. Everything that does not depend on k belongs here:
   !> the caller sets an instance up once and reuses it at every k-point, because the
   !> k dependence arrives with the coefficients and init_mat clears the result while
   !> keeping the allocation.
   !>
   !> A site is not summed over here. An instance covers one site, and whoever wants
   !> the total adds them up -- which is only a plain sum for an operator that is
   !> spin-diagonal with a frame-invariant trace. If yours is not, the caller has to
   !> rotate each site from its local frame first.
   SUBROUTINE melem_template_init(this, atoms, ntyp, iat, channel)
      CLASS(t_matelements_template), INTENT(INOUT) :: this
      TYPE(t_atoms), TARGET, INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: ntyp, iat
      INTEGER, INTENT(IN), OPTIONAL :: channel   !> present: serve this channel alone

      IF (ntyp < 1 .OR. ntyp > atoms%ntype) &
         CALL judft_bug("melem_template_init: the atom type is out of range")
      IF (iat < 1 .OR. iat > atoms%neq(ntyp)) &
         CALL judft_bug("melem_template_init: the equivalent atom is out of range")

      !> The coefficients always arrive for both spin components. What differs is what to
      !> do with them: over a spinor the two are one state and are summed, while two
      !> collinear channels are two states and each has its own value. Declaring it here
      !> lets the caller hand over whatever it has.
      this%channel = 0
      IF (PRESENT(channel)) THEN
         IF (channel < 1 .OR. channel > 2) CALL judft_bug("melem_template_init: the spin channel is 1 or 2")
         this%channel = channel
      END IF

      !> THE THREE VARIABLES THAT DECIDE THE SHAPE OF THE RESULT.
      !> The two spinor flags decide whether init_mat hands out a single block or a 2x2
      !> one in spin space: nsp = MERGE(2, 1, spinorwavefcts .OR. spinoroperator). An
      !> operator that does not act in spin space leaves both .FALSE. and gets one block.
      this%spinoroperator = .FALSE.
      this%spinorwavefcts = .FALSE.
      !> Cartesian components, declared here so that a caller which does not know which
      !> operator it is holding still allocates the right store. n_alpha = 1 means comp
      !> is unused and the result goes in mat instead. A component-carrying operator
      !> CANNOT be distributed over the eigenvector sub-communicator: comp has no
      !> distributed counterpart, and init_mat refuses the combination.
      this%n_alpha = 3

      this%atoms => atoms
      this%ntyp  = ntyp
      this%iat   = iat
   END SUBROUTINE melem_template_init

   SUBROUTINE melem_template_calc(this, zmat, abc, radfun, usdus)
      CLASS(t_matelements_template), INTENT(INOUT) :: this
      !> The state at this k-point in as few matrices as it takes: ONE when it is a whole
      !> spinor, TWO when the records are independent spin channels. SIZE(zmat) is
      !> therefore NOT jspins, and a consumer that addresses a spin block by row offset
      !> needs the first case. Unused here: a muffin-tin quantity works on abc alone.
      TYPE(t_mat),    INTENT(IN) :: zmat(:)
      TYPE(t_abc),    INTENT(IN) :: abc(:, :) !> (2 spin, ntype) local-frame coefficients
      TYPE(t_radfun), INTENT(IN) :: radfun(:) !> (ntype)
      TYPE(t_usdus),  INTENT(IN) :: usdus     !> unused: the radial integrals are in radfun

      INTEGER :: nb, i, j, l, ll1, mm, lm, n_r, n_r2, s, s_lo, s_hi, slot(2)
      REAL    :: w
      COMPLEX :: ovl, ang(3), acc(3)

      !> The guards are part of the contract, not decoration. Each one turns a shape
      !> mismatch into a message instead of into plausible numbers.
      IF (.NOT. ALLOCATED(this%mat)) &
         CALL judft_bug("melem_template_calc: the result matrix is not allocated")
      IF (SIZE(this%mat, 1) /= 1 .OR. SIZE(this%mat, 2) /= 1) &
         CALL judft_bug("melem_template_calc: this operator is spin-diagonal, so the result is a single block")
      IF (SIZE(abc, 1) /= 2 .OR. SIZE(abc, 2) /= this%atoms%ntype) &
         CALL judft_bug("melem_template_calc: the coefficients must have shape (2,ntype)")

      nb = SIZE(abc(1, 1)%cof, 1)
      IF (this%mat(1, 1)%matsize1 /= nb) &
         CALL judft_bug("melem_template_calc: the matrix size does not match the abc coefficients")

      !> One channel is one state, so only its own component is read. A spinor is one
      !> state made of two, so both are, and they are summed.
      IF (this%channel > 0) THEN
         s_lo = this%channel; s_hi = this%channel
      ELSE
         s_lo = 1; s_hi = 2
      END IF

      !> radfun%integral is allocated (.,.,.,jspins,jspins), so a component has to be
      !> clamped to the sets that exist. Resolved ONCE here: inside the loop it would be
      !> the same number fetched a million times.
      slot(1) = radial_slot(radfun, 1)
      slot(2) = radial_slot(radfun, 2)

      DO j = 1, nb                     ! ket band
         DO i = 1, nb                  ! bra band
            acc = CMPLX(0.0, 0.0)
            DO l = 0, this%atoms%lmax(this%ntyp)
               ll1 = l*(l + 1)
               DO mm = -l, l
                  lm = ll1 + mm
                  DO s = s_lo, s_hi    ! spin-diagonal: a spinor sums, a channel stands alone
                     DO n_r = 1, abc(s, this%ntyp)%n_r(l)
                        DO n_r2 = 1, abc(s, this%ntyp)%n_r(l)
                           w = radfun(this%ntyp)%integral(n_r, n_r2, l, slot(s), slot(s))
                           ovl = abc(s, this%ntyp)%cof(i, lm,  n_r,  this%iat) &
                                 * CONJG(abc(s, this%ntyp)%cof(j, lm, n_r2, this%iat)) * w

                           !> ------------------------------------------------------------
                           !> HERE GOES YOUR OPERATOR. What follows is the ANGULAR factor.
                           !> For the muffin-tin overlap it is 1 and does not depend on the
                           !> component, so the three come out equal -- which is exactly the
                           !> line to replace.
                           !>
                           !> In m_types_matelements_orbital the same place carries m for
                           !> L_z and the sqrt((l-+m)(l+-m+1)) of L_+ and L_-, and the
                           !> couplings there reach lm+1 and lm-1, so the ket index is NOT
                           !> the same lm as the bra. An operator that changes l or m reads
                           !> a NEIGHBOURING coefficient, and the bounds of that neighbour
                           !> have to be guarded (IF (mm < l), IF (mm > -l), ...).
                           !> ------------------------------------------------------------
                           ang(1) = CMPLX(1.0, 0.0)
                           ang(2) = CMPLX(1.0, 0.0)
                           ang(3) = CMPLX(1.0, 0.0)

                           acc(1:3) = acc(1:3) + ang(1:3)*ovl
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            !> Assigned and NOT accumulated: one instance covers one component of one
            !> atom, so each element gets exactly one contribution. Adding onto a cleared
            !> matrix would also turn a -0.0 result into +0.0, and the sign of zero is
            !> visible in the exported file.
            this%comp(i, j, 1:3) = acc(1:3)
         END DO
      END DO
   END SUBROUTINE melem_template_calc

END MODULE m_types_matelements_template
