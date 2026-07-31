!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Spin matrix elements <psi_m | sigma | psi_n> at one k-point, as the four 2x2
!>  spinor blocks. The three Pauli components follow from the blocks, so they are
!>  not stored:  sigma_x = M(1,2)+M(2,1),  sigma_y = -i(M(1,2)-M(2,1)),
!>  sigma_z = M(1,1)-M(2,2).
!>
!>  The operator has an interstitial and a muffin-tin part. The muffin-tin part is
!>  accumulated site by site, each site rotated from its own local frame to the
!>  global one before it is added: an antiferromagnet (beta = pi on one site) would
!>  otherwise add two local "+" moments and its total would not vanish.
MODULE m_types_matelements_spin
   USE m_types_matelements
   USE m_types_mat
   USE m_types_abc
   USE m_types_radfun
   USE m_types_usdus
   USE m_types_atoms
   USE m_types_stars
   USE m_types_lapw
   USE m_types_nococonv
   USE m_constants, ONLY: ImagUnit
   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_matelements) :: t_matelements_spin
      TYPE(t_atoms),    POINTER :: atoms    => NULL()
      TYPE(t_stars),    POINTER :: stars    => NULL()
      TYPE(t_lapw),     POINTER :: lapw     => NULL()
      TYPE(t_nococonv), POINTER :: nococonv => NULL()
   CONTAINS
      PROCEDURE :: init                 => init
      PROCEDURE :: calc_matrix_elements => calc_matrix_elements
   END TYPE t_matelements_spin

   PUBLIC :: t_matelements_spin

CONTAINS

   SUBROUTINE init(this, atoms, stars, lapw, nococonv)
      CLASS(t_matelements_spin), INTENT(INOUT) :: this
      TYPE(t_atoms),    TARGET,  INTENT(IN)    :: atoms
      TYPE(t_stars),    TARGET,  INTENT(IN)    :: stars
      TYPE(t_lapw),     TARGET,  INTENT(IN)    :: lapw     !> changes with k
      TYPE(t_nococonv), TARGET,  INTENT(IN)    :: nococonv

      !> sigma couples the two spin channels, so the result is a 2x2 block matrix
      this%spinoroperator = .TRUE.
      this%spinorwavefcts = .TRUE.

      this%atoms    => atoms
      this%stars    => stars
      this%lapw     => lapw
      this%nococonv => nococonv
   END SUBROUTINE init

   SUBROUTINE calc_matrix_elements(this, zmat, abc, radfun, usdus)
      CLASS(t_matelements_spin), INTENT(INOUT) :: this
      TYPE(t_mat),    INTENT(IN) :: zmat(:)   !> the 2N spinor in one matrix
      TYPE(t_abc),    INTENT(IN) :: abc(:,:)  !> (2,ntype)
      TYPE(t_radfun), INTENT(IN) :: radfun(:) !> (ntype)
      TYPE(t_usdus),  INTENT(IN) :: usdus     !> unused, the radial integrals are in radfun

      COMPLEX, ALLOCATABLE :: oi(:,:,:,:)     ! (nb,nb,2,2) interstitial spin blocks
      COMPLEX :: loc(2,2), glo(2,2), cx, cy, cz, gx, gy, gz, trc
      REAL    :: ca, sa, cb, sb
      INTEGER :: nb, io_dn, i, j, ntyp, iat, l, ll1, mm, lm, n_r, n_r2
      INTEGER :: js2, i1, j1
      LOGICAL :: l_rot

      IF (.NOT.ALLOCATED(this%mat)) &
         CALL judft_bug("calc_matrix_elements: the result matrix is not allocated")
      IF (SIZE(this%mat,1) /= 2 .OR. SIZE(this%mat,2) /= 2) &
         CALL judft_bug("calc_matrix_elements: the matrix must be a 2x2 spinor matrix")
      IF (SIZE(abc,1) /= 2 .OR. SIZE(abc,2) /= this%atoms%ntype) &
         CALL judft_bug("calc_matrix_elements: the abc coefficients must have shape (2,ntype)")

      nb = SIZE(abc(1,1)%cof, 1)
      IF (this%mat(1,1)%matsize1 /= nb) &
         CALL judft_bug("calc_matrix_elements: the matrix size does not match the abc coefficients")

      !> The radial set is only input%jspins wide, so a single set serves both
      !> spinor components.
      js2 = MERGE(1, 2, SIZE(radfun(1)%integral, 4) < 2)

      ! ---- interstitial: the spin-down rows of the spinor start at io_dn ----
      io_dn = this%lapw%nv(1) + this%atoms%nlotot
      ALLOCATE(oi(nb, nb, 2, 2))
      CALL overlap_interstitial(this%stars, this%lapw, zmat(1), 0,     0,     oi(:,:,1,1))
      CALL overlap_interstitial(this%stars, this%lapw, zmat(1), io_dn, io_dn, oi(:,:,2,2))
      CALL overlap_interstitial(this%stars, this%lapw, zmat(1), 0,     io_dn, oi(:,:,1,2))
      CALL overlap_interstitial(this%stars, this%lapw, zmat(1), io_dn, 0,     oi(:,:,2,1))
      DO j1 = 1, 2
         DO i1 = 1, 2
            this%mat(i1,j1)%data_c(:,:) = this%mat(i1,j1)%data_c(:,:) + oi(:,:,i1,j1)
         END DO
      END DO
      DEALLOCATE(oi)

      ! ---- muffin-tin, site by site, each rotated into the global frame ----
      DO j = 1, nb
         DO i = 1, nb
            DO ntyp = 1, this%atoms%ntype
               !> A site whose frame is already the global one needs no rotation, and
               !> going through the Pauli components and back would only cost precision.
               l_rot = ABS(this%nococonv%alph(ntyp)) > 1.0e-14 &
                  .OR. ABS(this%nococonv%beta(ntyp)) > 1.0e-14
               ca = COS(this%nococonv%alph(ntyp)); sa = SIN(this%nococonv%alph(ntyp))
               cb = COS(this%nococonv%beta(ntyp)); sb = SIN(this%nococonv%beta(ntyp))
               DO iat = 1, this%atoms%neq(ntyp)
                  loc = CMPLX(0.0, 0.0)
                  DO l = 0, this%atoms%lmax(ntyp)
                     ll1 = l*(l + 1)
                     DO mm = -l, l
                        lm = ll1 + mm
                        DO n_r = 1, abc(1,ntyp)%n_r(l)
                           DO n_r2 = 1, abc(1,ntyp)%n_r(l)
                              loc(1,1) = loc(1,1) + abc(1,ntyp)%cof(i,lm,n_r,iat) &
                                       * CONJG(abc(1,ntyp)%cof(j,lm,n_r2,iat)) &
                                       * radfun(ntyp)%integral(n_r,n_r2,l,1,1)
                              loc(2,2) = loc(2,2) + abc(2,ntyp)%cof(i,lm,n_r,iat) &
                                       * CONJG(abc(2,ntyp)%cof(j,lm,n_r2,iat)) &
                                       * radfun(ntyp)%integral(n_r,n_r2,l,js2,js2)
                              loc(1,2) = loc(1,2) + abc(1,ntyp)%cof(i,lm,n_r,iat) &
                                       * CONJG(abc(2,ntyp)%cof(j,lm,n_r2,iat)) &
                                       * radfun(ntyp)%integral(n_r,n_r2,l,1,js2)
                              loc(2,1) = loc(2,1) + abc(2,ntyp)%cof(i,lm,n_r,iat) &
                                       * CONJG(abc(1,ntyp)%cof(j,lm,n_r2,iat)) &
                                       * radfun(ntyp)%integral(n_r,n_r2,l,js2,1)
                           END DO
                        END DO
                     END DO
                  END DO
                  IF (l_rot) THEN
                     ! local Pauli components, rotated local->global by R_z(alpha) R_y(beta),
                     ! then back to a block; the trace is invariant under the rotation
                     cx = loc(1,2) + loc(2,1)
                     cy = -ImagUnit * (loc(1,2) - loc(2,1))
                     cz = loc(1,1) - loc(2,2)
                     gx =  ca*cb*cx - sa*cy + ca*sb*cz
                     gy =  sa*cb*cx + ca*cy + sa*sb*cz
                     gz = -sb*cx           + cb*cz
                     trc = loc(1,1) + loc(2,2)
                     glo(1,1) = 0.5 * (trc + gz)
                     glo(2,2) = 0.5 * (trc - gz)
                     glo(1,2) = 0.5 * (gx + ImagUnit*gy)
                     glo(2,1) = 0.5 * (gx - ImagUnit*gy)
                  ELSE
                     glo = loc
                  END IF
                  DO j1 = 1, 2
                     DO i1 = 1, 2
                        this%mat(i1,j1)%data_c(i,j) = this%mat(i1,j1)%data_c(i,j) + glo(i1,j1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE calc_matrix_elements

   !> Interstitial overlap of two spinor components of one eigenvector at one k,
   !>    ovl(m,n) = <psi^a_m | psi^b_n>_interstitial ,
   !> with ioff_a and ioff_b selecting the rows of zmat that carry each component.
   !> The step function supplies the overlap of two plane waves over the interstitial,
   !> Theta(G_n - G_m); ig maps that difference to its star and rgphs its phase.
   SUBROUTINE overlap_interstitial(stars, lapw, zmat, ioff_a, ioff_b, ovl)
      TYPE(t_stars), INTENT(IN)  :: stars
      TYPE(t_lapw),  INTENT(IN)  :: lapw
      TYPE(t_mat),   INTENT(IN)  :: zmat
      INTEGER,       INTENT(IN)  :: ioff_a, ioff_b
      COMPLEX,       INTENT(OUT) :: ovl(:,:)

      COMPLEX, ALLOCATABLE :: stepf(:,:), zstep(:,:)
      INTEGER :: nv, nb, i, j, i1, i2, i3, in

      IF (zmat%l_real) &
         CALL judft_bug("overlap_interstitial: a spinor cannot be held in a real matrix")

      nv = lapw%nv(1)
      nb = SIZE(ovl, 1)
      ALLOCATE(stepf(nv, nv), source=CMPLX(0.0, 0.0))
      ALLOCATE(zstep(nv, nb))

      DO i = 1, nv
         DO j = 1, nv
            i1 = lapw%k1(j,1) - lapw%k1(i,1)
            i2 = lapw%k2(j,1) - lapw%k2(i,1)
            i3 = lapw%k3(j,1) - lapw%k3(i,1)
            in = stars%ig(i1, i2, i3)
            IF (in == 0) CYCLE
            stepf(j,i) = CONJG(stars%rgphs(i1,i2,i3) * stars%ustep(in))
         END DO
      END DO

      CALL zgemm('T', 'N', nv, nb, nv, CMPLX(1.0,0.0), stepf, nv, &
                 zmat%data_c(1+ioff_b, 1), zmat%matsize1, CMPLX(0.0,0.0), zstep, nv)
      zstep = CONJG(zstep)
      CALL zgemm('T', 'N', nb, nb, nv, CMPLX(1.0,0.0), zmat%data_c(1+ioff_a, 1), zmat%matsize1, &
                 zstep, nv, CMPLX(0.0,0.0), ovl, nb)
   END SUBROUTINE overlap_interstitial

END MODULE m_types_matelements_spin
