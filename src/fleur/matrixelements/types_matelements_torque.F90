!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Exchange-torque matrix elements <psi_m | B_xc | psi_n> at one k-point, as the
!>  four 2x2 spinor blocks, with B_xc = (V_up - V_dn)/2 the exchange field.
!>
!>  It is the same object the classic route writes as updown.perpmag, and the
!>  torque of Freimuth, Mokrousov et al. follows from the blocks the way the Pauli
!>  components do: with the moment along the local z, T = -mu_B sigma x Omega_xc
!>  has its two non-vanishing components in the spin off-diagonal blocks, so
!>  melem_pauli_from_blocks applies here unchanged and nothing extra is stored.
!>
!>  Two things set it apart from the spin operator it is otherwise a copy of:
!>
!>    - the weight. The spin operator contracts radfun%integral, a plain overlap,
!>      which the orthogonality of the Y_lm makes diagonal in lm. B_xc is not
!>      spherical, so the weight is t_melem_beff, which couples lm to lm'.
!>
!>    - the regions. B_xc is a muffin-tin quantity here: the interstitial and the
!>      vacuum are not summed, which is what the original does too. That makes the
!>      operator meaningless on a film, and init says so rather than returning a
!>      number that silently misses a region.
MODULE m_types_matelements_torque
   USE m_types_matelements
   USE m_types_mat
   USE m_types_abc
   USE m_types_radfun
   USE m_types_melem_beff, ONLY: t_melem_beff
   USE m_types_input
   USE m_types_atoms
   USE m_types_sym
   USE m_types_sphhar
   USE m_types_potden
   USE m_types_nococonv
   USE m_constants, ONLY: ImagUnit
   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_matelements) :: t_matelements_torque
      TYPE(t_atoms),    POINTER :: atoms    => NULL()
      TYPE(t_nococonv), POINTER :: nococonv => NULL()
      TYPE(t_sym),      POINTER :: sym      => NULL()
      TYPE(t_sphhar),   POINTER :: sphhar   => NULL()
      TYPE(t_potden),   POINTER :: vtot     => NULL()
      !> The weight. It needs the radial functions, which only arrive with the states,
      !> so it is built on the first call rather than at init and kept afterwards: it
      !> does not depend on k.
      TYPE(t_melem_beff)        :: beff
   CONTAINS
      PROCEDURE :: init                 => torque_init
      PROCEDURE :: calc_matrix_elements => torque_calc_matrix_elements
   END TYPE t_matelements_torque

   PUBLIC :: t_matelements_torque

CONTAINS

   !> Everything this operator needs is in abc and in the potential: it never reads
   !> zmat, so it takes neither the basis nor the layout that says where in zmat a
   !> spin component lives.
   SUBROUTINE torque_init(this, atoms, sym, sphhar, nococonv, input, vtot)
      CLASS(t_matelements_torque), INTENT(INOUT) :: this
      TYPE(t_atoms),    TARGET,  INTENT(IN) :: atoms
      TYPE(t_sym),      TARGET,  INTENT(IN) :: sym
      TYPE(t_sphhar),   TARGET,  INTENT(IN) :: sphhar
      TYPE(t_nococonv), TARGET,  INTENT(IN) :: nococonv
      TYPE(t_input),             INTENT(IN) :: input
      TYPE(t_potden),   TARGET,  INTENT(IN) :: vtot

      IF (input%film) &
         CALL judft_error('the exchange torque has no vacuum contribution and cannot be '// &
                          'asked for on a film', &
                          hint='it sums the muffin tins only; orbital and spinOrbit are '// &
                               'muffin-tin quantities too and do work on a film', &
                          calledby='types_matelements_torque%init')

      !> B_xc couples the two spin channels, so the result is a 2x2 block matrix
      this%spinoroperator = .TRUE.
      this%spinorwavefcts = .TRUE.

      this%atoms    => atoms
      this%nococonv => nococonv
      this%sym      => sym
      this%sphhar   => sphhar
      this%vtot     => vtot
   END SUBROUTINE torque_init

   SUBROUTINE torque_calc_matrix_elements(this, zmat, abc, radfun)
      CLASS(t_matelements_torque), INTENT(INOUT) :: this
      TYPE(t_mat),    INTENT(IN) :: zmat(:)  !> unused: the deferred interface carries it
      TYPE(t_abc),    INTENT(IN) :: abc(:,:)  !> (2,ntype)
      TYPE(t_radfun), INTENT(IN) :: radfun(:) !> (ntype)

      COMPLEX, ALLOCATABLE :: blk(:,:,:,:)   ! (nb,nb,2,2) one site, local frame
      COMPLEX, ALLOCATABLE :: w(:,:,:)       ! (nb,0:nlm-1,nrr) the half-contracted weight
      COMPLEX, ALLOCATABLE :: cx(:,:), cy(:,:), cz(:,:), trc(:,:)
      COMPLEX, ALLOCATABLE :: gx(:,:), gy(:,:), gz(:,:)
      REAL    :: ca, sa_, cb, sb_
      INTEGER :: nb, ntyp, iat, nlm, nrr, lma, lmb, ia, ib, sa, sb, i1, j1
      LOGICAL :: l_rot

      IF (.NOT.ALLOCATED(this%mat)) &
         CALL judft_bug("calc_matrix_elements: the result matrix is not allocated")
      IF (SIZE(this%mat,1) /= 2 .OR. SIZE(this%mat,2) /= 2) &
         CALL judft_bug("calc_matrix_elements: the matrix must be a 2x2 spinor matrix")
      IF (SIZE(abc,1) /= 2 .OR. SIZE(abc,2) /= this%atoms%ntype) &
         CALL judft_bug("calc_matrix_elements: the abc coefficients must have shape (2,ntype)")

      !> The weight does not depend on k, so it is built once, here, where radfun is
      !> in scope for the first time.
      IF (.NOT.ALLOCATED(this%beff%g)) &
         CALL this%beff%init(this%atoms, this%sym, this%sphhar, this%vtot, radfun)

      nb = SIZE(abc(1,1)%cof, 1)
      IF (this%mat(1,1)%matsize1 /= nb) &
         CALL judft_bug("calc_matrix_elements: the matrix size does not match the abc coefficients")

      ALLOCATE(blk(nb, nb, 2, 2))
      ALLOCATE(cx(nb,nb), cy(nb,nb), cz(nb,nb), trc(nb,nb))
      ALLOCATE(gx(nb,nb), gy(nb,nb), gz(nb,nb))

      ! ---- muffin-tin, site by site, each rotated into the global frame ----
      DO ntyp = 1, this%atoms%ntype
         nlm = SIZE(abc(1,ntyp)%cof, 2)
         nrr = SIZE(abc(1,ntyp)%cof, 3)
         IF (nlm - 1 > this%beff%lmd .OR. nrr > this%beff%n_rd) &
            CALL judft_bug("calc_matrix_elements: the weight table is smaller than the abc &
                           &coefficients it has to contract")

         !> A site whose frame is already the global one needs no rotation, and going
         !> through the Pauli components and back would only cost precision.
         !>
         !> The azimuth enters with a MINUS sign, and that sign and the single
         !> conjugation at the end of this routine belong together, exactly as in the
         !> spin operator this is a copy of:
         !>     conj( R_z(-alpha) [ conj(B) ] ) = R_z(+alpha) [ B ] .
         !> Changing either one alone breaks it.
         l_rot = ABS(this%nococonv%alph(ntyp)) > 1.0e-14 &
            .OR. ABS(this%nococonv%beta(ntyp)) > 1.0e-14
         ca = COS(this%nococonv%alph(ntyp)); sa_ = -SIN(this%nococonv%alph(ntyp))
         cb = COS(this%nococonv%beta(ntyp)); sb_ =  SIN(this%nococonv%beta(ntyp))

         ALLOCATE(w(nb, 0:nlm-1, nrr))
         DO iat = 1, this%atoms%neq(ntyp)
            blk = CMPLX(0.0, 0.0)

            !> Contract in two passes over the composite index (lm, slot). Done band by
            !> band the two lm sums would cost nb^2 * nlm^2 * nrr^2; this way the band
            !> index rides along and it costs nb*nlm^2*nrr^2 + nb^2*nlm*nrr.
            DO sb = 1, 2
               DO sa = 1, 2
                  w = CMPLX(0.0, 0.0)
                  DO ib = 1, nrr
                     DO lmb = 0, nlm-1
                        DO ia = 1, nrr
                           DO lma = 0, nlm-1
                              w(:, lmb, ib) = w(:, lmb, ib)                    &
                                 + abc(sa,ntyp)%cof(:, lma, ia, iat)           &
                                 * CONJG(this%beff%g(lma, ia, lmb, ib, sa, sb, ntyp))
                           END DO
                        END DO
                     END DO
                  END DO
                  DO ib = 1, nrr
                     DO lmb = 0, nlm-1
                        blk(:,:,sa,sb) = blk(:,:,sa,sb)                        &
                           + SPREAD(w(:, lmb, ib), 2, nb)                      &
                           * SPREAD(CONJG(abc(sb,ntyp)%cof(:, lmb, ib, iat)), 1, nb)
                     END DO
                  END DO
               END DO
            END DO

            IF (l_rot) THEN
               ! local Pauli components, rotated local->global by R_z(-alpha) R_y(beta),
               ! then back to blocks; the trace is invariant under the rotation
               cx = blk(:,:,1,2) + blk(:,:,2,1)
               cy = -ImagUnit * (blk(:,:,1,2) - blk(:,:,2,1))
               cz = blk(:,:,1,1) - blk(:,:,2,2)
               gx =  ca*cb*cx - sa_*cy + ca*sb_*cz
               gy =  sa_*cb*cx + ca*cy + sa_*sb_*cz
               gz = -sb_*cx            + cb*cz
               trc = blk(:,:,1,1) + blk(:,:,2,2)
               blk(:,:,1,1) = 0.5 * (trc + gz)
               blk(:,:,2,2) = 0.5 * (trc - gz)
               blk(:,:,1,2) = 0.5 * (gx + ImagUnit*gy)
               blk(:,:,2,1) = 0.5 * (gx - ImagUnit*gy)
            END IF

            DO j1 = 1, 2
               DO i1 = 1, 2
                  this%mat(i1,j1)%data_c(:,:) = this%mat(i1,j1)%data_c(:,:) + blk(:,:,i1,j1)
               END DO
            END DO
         END DO
         DEALLOCATE(w)
      END DO

      !> The loop above forms c_i * conj(c_j) against conj(g), the conjugate of the
      !> matrix element <i|B_xc|j> this type promises, so one conjugation of the
      !> assembled blocks restores the promised convention. It is the same convention
      !> the spin operator uses, and the azimuthal sign above depends on it.
      DO j1 = 1, 2
         DO i1 = 1, 2
            this%mat(i1,j1)%data_c = CONJG(this%mat(i1,j1)%data_c)
         END DO
      END DO

      DEALLOCATE(blk, cx, cy, cz, trc, gx, gy, gz)
   END SUBROUTINE torque_calc_matrix_elements

END MODULE m_types_matelements_torque
