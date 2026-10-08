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
   USE m_types_spinor_layout, ONLY: t_spinor_layout, radial_slot, LAYOUT_SPINOR, LAYOUT_CHANNELS
   USE m_types_input
   USE m_types_noco
   USE m_types_usdus
   USE m_types_atoms
   USE m_types_stars
   USE m_types_lapw
   USE m_types_nococonv
   USE m_melem_overlap, ONLY: melem_overlap_interstitial
  USE m_types_melem_vacabc, ONLY: t_melem_vacabc
  USE m_melem_mmkb_vac, ONLY: melem_mmkb_vac
  USE m_types_vacuum
  USE m_types_cell
  USE m_types_enpara
  USE m_types_potden
   USE m_constants, ONLY: ImagUnit
   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_matelements) :: t_matelements_spin
      TYPE(t_atoms),    POINTER :: atoms    => NULL()
      TYPE(t_stars),    POINTER :: stars    => NULL()
      TYPE(t_lapw),     POINTER :: lapw     => NULL()
      TYPE(t_nococonv), POINTER :: nococonv => NULL()
      !> Only a film has a third region. They are pointers and not copies for the same
      !> reason the four above are, and l_film is what says they were handed over: a bulk
      !> caller leaves them out and the vacuum half is then simply not there.
      TYPE(t_vacuum),   POINTER :: vacuum   => NULL()
      TYPE(t_cell),     POINTER :: cell     => NULL()
      TYPE(t_enpara),   POINTER :: enpara   => NULL()
      TYPE(t_potden),   POINTER :: vtot     => NULL()
      LOGICAL                   :: l_film   = .FALSE.
      TYPE(t_spinor_layout)     :: layout
   CONTAINS
      PROCEDURE :: init                 => init
      PROCEDURE :: calc_matrix_elements => calc_matrix_elements
   END TYPE t_matelements_spin

   PUBLIC :: t_matelements_spin

CONTAINS

   SUBROUTINE init(this, atoms, stars, lapw, nococonv, input, noco, &
                   vacuum, cell, enpara, vtot)
      CLASS(t_matelements_spin), INTENT(INOUT) :: this
      TYPE(t_atoms),    TARGET,  INTENT(IN)    :: atoms
      TYPE(t_stars),    TARGET,  INTENT(IN)    :: stars
      TYPE(t_lapw),     TARGET,  INTENT(IN)    :: lapw     !> changes with k
      TYPE(t_nococonv), TARGET,  INTENT(IN)    :: nococonv
      TYPE(t_input),             INTENT(IN)    :: input
      TYPE(t_noco),             INTENT(IN)     :: noco
      !> What the vacuum expansion needs. Optional so that a bulk caller says nothing and
      !> gets the two regions it has, rather than having to pass four unused objects.
      TYPE(t_vacuum), TARGET, OPTIONAL, INTENT(IN) :: vacuum
      TYPE(t_cell),   TARGET, OPTIONAL, INTENT(IN) :: cell
      TYPE(t_enpara), TARGET, OPTIONAL, INTENT(IN) :: enpara
      TYPE(t_potden), TARGET, OPTIONAL, INTENT(IN) :: vtot

      !> sigma couples the two spin channels, so the result is a 2x2 block matrix
      this%spinoroperator = .TRUE.
      this%spinorwavefcts = .TRUE.

      this%atoms    => atoms
      this%stars    => stars
      this%lapw     => lapw
      this%nococonv => nococonv

      !> All four or none: a film that is missing one of them would silently lose the
      !> vacuum half, which is exactly the failure this is meant to remove.
      this%l_film = .FALSE.
      IF (input%film) THEN
         IF (.NOT.(PRESENT(vacuum) .AND. PRESENT(cell) .AND. PRESENT(enpara) .AND. PRESENT(vtot))) &
            CALL judft_error('t_matelements_spin: on a film the spin operator needs the vacuum', &
                             hint='pass vacuum, cell, enpara and vtot to init; without them the '// &
                                  'four blocks would hold the muffin tins and the interstitial only', &
                             calledby='types_matelements_spin%init')
         this%vacuum => vacuum; this%cell => cell
         this%enpara => enpara; this%vtot => vtot
         this%l_film = .TRUE.
      END IF

      !> The caller hands over a whole spinor either way: one 2N record when the
      !> Hamiltonian was non-collinear, two records stacked by the caller otherwise.
      !> Saying so is what makes row_dn meaningful here.
      CALL this%layout%init(input, noco, lapw, atoms, &
                            l_both_spinors=(noco%l_soc .AND. .NOT.noco%l_noco))
   END SUBROUTINE init

   SUBROUTINE calc_matrix_elements(this, zmat, abc, radfun, usdus)
      CLASS(t_matelements_spin), INTENT(INOUT) :: this
      !> The state at this k, either as one 2N matrix holding the whole spinor or as the
      !> two channels of a collinear calculation, one matrix each. Which of the two is read
      !> off its size, and it decides where the spin-down rows are: inside the single matrix
      !> at row_dn, or at the top of the second one.
      TYPE(t_mat),    INTENT(IN) :: zmat(:)
      TYPE(t_abc),    INTENT(IN) :: abc(:,:)  !> (2,ntype)
      TYPE(t_radfun), INTENT(IN) :: radfun(:) !> (ntype)
      TYPE(t_usdus),  INTENT(IN) :: usdus     !> unused, the radial integrals are in radfun

      COMPLEX, ALLOCATABLE :: oi(:,:,:,:)     ! (nb,nb,2,2) interstitial spin blocks
      COMPLEX, ALLOCATABLE :: ov(:,:)         ! (nb,nb) one vacuum block at a time
      TYPE(t_melem_vacabc) :: vacs(2)         ! the expansion of each spin channel
      COMPLEX :: loc(2,2), glo(2,2), cx, cy, cz, gx, gy, gz, trc
      REAL    :: ca, sa, cb, sb
      INTEGER :: nb, i, j, ntyp, iat, l, ll1, mm, lm, n_r, n_r2
      INTEGER :: js2, i1, j1
      INTEGER :: iu, id, off_u, off_d
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

      js2 = radial_slot(radfun, 2)

      !> Where each spin component is to be found. Stacked, both live in one matrix and the
      !> down rows begin at row_dn; as two channels, each has its own matrix and starts at
      !> the top. The rest of this routine does not care which of the two it was.
      SELECT CASE (SIZE(zmat))
      CASE (1)
         IF (this%layout%layout /= LAYOUT_SPINOR) &
            CALL judft_bug("calc_matrix_elements: one matrix holds a whole spinor, and the &
                           &layout of this k-point does not say so")
         iu = 1; id = 1; off_u = 0; off_d = this%layout%row_dn
      CASE (2)
         IF (this%layout%layout /= LAYOUT_CHANNELS) &
            CALL judft_bug("calc_matrix_elements: two matrices are two spin channels, and the &
                           &layout of this k-point does not say so")
         iu = 1; id = 2; off_u = 0; off_d = 0
      CASE DEFAULT
         CALL judft_bug("calc_matrix_elements: the state arrives either as one spinor or as &
                        &two channels")
      END SELECT

      ! ---- interstitial, one call per spin block ----
      ALLOCATE(oi(nb, nb, 2, 2))
      CALL melem_overlap_interstitial(this%stars, this%lapw, this%lapw, zmat(iu), zmat(iu), &
                                      off_u, off_u, oi(:,:,1,1))
      CALL melem_overlap_interstitial(this%stars, this%lapw, this%lapw, zmat(id), zmat(id), &
                                      off_d, off_d, oi(:,:,2,2))
      CALL melem_overlap_interstitial(this%stars, this%lapw, this%lapw, zmat(iu), zmat(id), &
                                      off_u, off_d, oi(:,:,1,2))
      CALL melem_overlap_interstitial(this%stars, this%lapw, this%lapw, zmat(id), zmat(iu), &
                                      off_d, off_u, oi(:,:,2,1))
      DO j1 = 1, 2
         DO i1 = 1, 2
            this%mat(i1,j1)%data_c(:,:) = this%mat(i1,j1)%data_c(:,:) + oi(:,:,i1,j1)
         END DO
      END DO
      DEALLOCATE(oi)

      ! ---- vacuum, the third region a film has ----
      !> One expansion per spin channel, and a block is the bra channel against the ket
      !> channel: o_ab = <phi^a|phi^b>_vac, the same pairing the interstitial above uses.
      !> gb is zero because both sides are the same k.
      !>
      !> It accumulates in the same c_i conj(c_j) convention as the two halves above, so it
      !> goes in BEFORE the single conjugation at the end of this routine, not after.
      IF (this%l_film) THEN
         ALLOCATE(ov(nb, nb))
         CALL vacs(1)%calc(this%vacuum, this%cell, this%enpara, this%vtot, this%lapw, &
                           radial_slot(radfun, 1), zmat(iu), nb, ioff=off_u)
         CALL vacs(2)%calc(this%vacuum, this%cell, this%enpara, this%vtot, this%lapw, &
                           radial_slot(radfun, 2), zmat(id), nb, ioff=off_d)
         DO j1 = 1, 2
            DO i1 = 1, 2
               ov = CMPLX(0.0, 0.0)
               CALL melem_mmkb_vac(vacs(i1), vacs(j1), [0, 0, 0], ov)
               this%mat(i1,j1)%data_c(:,:) = this%mat(i1,j1)%data_c(:,:) + ov(:,:)
            END DO
         END DO
         DEALLOCATE(ov)
      END IF

      ! ---- muffin-tin, site by site, each rotated into the global frame ----
      DO j = 1, nb
         DO i = 1, nb
            DO ntyp = 1, this%atoms%ntype
               !> A site whose frame is already the global one needs no rotation, and
               !> going through the Pauli components and back would only cost precision.
               !>
               !> The azimuth enters with a MINUS sign. A non-collinear basis already
               !> carries it in the exp(-i alpha/2), exp(+i alpha/2) phases of its two
               !> spinor components, so the cross-spin block o_ud arrives with the
               !> azimuthal phase already in it, and R_z(+alpha) on top would count that
               !> rotation twice instead of removing it. The polar rotation is unaffected.
               l_rot = ABS(this%nococonv%alph(ntyp)) > 1.0e-14 &
                  .OR. ABS(this%nococonv%beta(ntyp)) > 1.0e-14
               ca = COS(this%nococonv%alph(ntyp)); sa = -SIN(this%nococonv%alph(ntyp))
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
                     ! R_z(-alpha) R_y(beta).
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

      !> Both halves above form c_i * conj(c_j), the conjugate of the matrix element
      !> <i|O|j> this type promises, and they agree with each other, so one conjugation of
      !> the assembled blocks restores the promised convention.
      !>
      !> This conjugation and the sign of the azimuthal rotation belong together:
      !> conjugation flips the y component, and flipping y turns a rotation about z into
      !> its inverse, so
      !>     conj( R_z(-alpha) [ conj(B) ] ) = R_z(+alpha) [ B ] ,
      !> the canonical rotation. Changing either sign alone breaks it.
      DO j1 = 1, 2
         DO i1 = 1, 2
            this%mat(i1,j1)%data_c = CONJG(this%mat(i1,j1)%data_c)
         END DO
      END DO
   END SUBROUTINE calc_matrix_elements

END MODULE m_types_matelements_spin
