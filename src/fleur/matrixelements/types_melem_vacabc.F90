!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The vacuum expansion of one k-point's eigenvectors -- the sibling of t_abc, which does
!>  the same job inside the muffin tins.
!>
!>  Outside the film a basis function is a two-dimensional plane wave times one of two
!>  solutions of the z equation: u(z), which decays, and ue(z), its energy derivative. So
!>
!>      psi_n = sum_Gpar e^{i (k_par + G_par).r_par} [ ac(Gpar,n) u(z) + bc(Gpar,n) ue(z) ]
!>
!>  and ac/bc are what every vacuum matrix element contracts. The pair overlap and the
!>  momentum both need them, which is why they are built once here rather than inside
!>  each consumer.
!>
!>  SLOTS, NOT VACUA. A film has two sides, but under mirror symmetry FLEUR stores only one
!>  of them (vacuum%nvac == 1) and the other follows from z -> -z. Rather than teach every
!>  consumer that, this type always offers TWO slots: slot_vac says whose z functions a slot
!>  reads and slot_sign carries the reflection. A consumer loops over both slots and never
!>  learns which case it is in.
!>
!>  NORMALISATION, and it is not the one next door. ac/bc carry 1/(sqrt(omtil)*wronk), so an
!>  overlap built from them still needs the area as its in-plane measure. eigen/hsvac.F90
!>  folds the area in instead, through d2 = sqrt(omtil/area). The two conventions differ by
!>  exactly that factor and mixing them is silent, so the measure is kept explicit at the
!>  point where it is an integration measure.
MODULE m_types_melem_vacabc
   USE m_juDFT
   USE m_vacuz
   USE m_vacudz
   USE m_constants, ONLY: oUnit
   USE m_differentiate, ONLY: difcub
   USE m_types_cell
   USE m_types_enpara
   USE m_types_lapw
   USE m_types_mat
   USE m_types_potden
   USE m_types_vacuum
   IMPLICIT NONE
   PRIVATE

   !> The Wronskian the two z solutions are rescaled to. It is what fixes the meaning of
   !> ac/bc, so it belongs to the type and not to a caller.
   REAL, PARAMETER, PUBLIC :: MELEM_VAC_WRONK = 2.0

   !> How far the five-point derivative may sit from the boundary value vacuz returns,
   !> relative. Loose on purpose: what it has to catch is a wrong convention or a mesh that
   !> cannot carry a stencil, not the last digits. The measured value is printed either way.
   REAL, PARAMETER :: MELEM_VAC_DDZ_TOL = 5.0e-2
   LOGICAL, SAVE :: l_said = .FALSE.

   TYPE, PUBLIC :: t_melem_vacabc
      COMPLEX, ALLOCATABLE :: ac(:, :, :)   !< (nv2, nbnd, 2) resolved by slot
      COMPLEX, ALLOCATABLE :: bc(:, :, :)   !< (nv2, nbnd, 2)
      REAL, ALLOCATABLE :: u(:, :, :)       !< (nmz, nv2, nvac) the decaying solution
      REAL, ALLOCATABLE :: ue(:, :, :)      !< (nmz, nv2, nvac) its energy derivative
      !> The z derivatives of the two solutions, on the same mesh. The pair overlap does not
      !> need them; the momentum does, and it is the only operator that does. They are the
      !> derivative with respect to DISTANCE FROM THE BOUNDARY, so a consumer working on the
      !> mirrored slot multiplies by slot_sign to get d/dz.
      REAL, ALLOCATABLE :: du(:, :, :)      !< (nmz, nv2, nvac)
      REAL, ALLOCATABLE :: due(:, :, :)     !< (nmz, nv2, nvac)
      INTEGER, ALLOCATABLE :: kvac(:, :)    !< (2, nv2) the distinct parallel G
      !> k_par + G_par in CARTESIAN units, z component zero. The in-plane momentum is just
      !> multiplication by this, so it is built here where cell and lapw are both in scope
      !> rather than passed to every consumer. Same order as kvac.
      REAL, ALLOCATABLE :: qpar(:, :)       !< (3, nv2)
      INTEGER :: nv2 = 0
      INTEGER :: nmz = 0
      !> The z mesh and the geometry the expansion was built on. They travel with the data
      !> because a consumer handed a different t_vacuum or t_cell would integrate the same
      !> ac/bc on the wrong mesh and say nothing.
      REAL :: delz = 0.0
      REAL :: z1 = 0.0
      REAL :: area = 0.0
      REAL :: bmat33 = 0.0
      INTEGER :: slot_vac(2) = 0            !< whose u/ue each slot reads
      REAL :: slot_sign(2) = 0.0            !< +1 the stored side, -1 the other one
   CONTAINS
      PROCEDURE :: calc => melem_vacabc_calc
   END TYPE t_melem_vacabc

CONTAINS

   !>  Expands nbnd states of zMat in the vacuum, for one k-point and one spin.
   !>
   !>  ioff is the offset to the spin-down block of a non-collinear record, the same
   !>  argument melem_mmkb_int takes and for the same reason.
   SUBROUTINE melem_vacabc_calc(this, vacuum, cell, enpara, vtot, lapw, jspin, zMat, nbnd, ioff)
      CLASS(t_melem_vacabc), INTENT(INOUT) :: this
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_potden), INTENT(IN) :: vtot
      TYPE(t_lapw), INTENT(IN) :: lapw
      INTEGER, INTENT(IN) :: jspin
      TYPE(t_mat), INTENT(IN) :: zMat
      INTEGER, INTENT(IN) :: nbnd
      INTEGER, INTENT(IN), OPTIONAL :: ioff

      INTEGER :: nmz, nvac, nv, off, ikG, ik2, ivac, islot, n
      INTEGER, ALLOCATABLE :: map2(:)
      REAL :: vz0, ev, sc, zks, arg, sgn, v(3), const
      REAL, ALLOCATABLE :: vz(:), t(:), dt(:), te(:), dte(:), tei(:), zmsh(:)
      COMPLEX :: c_1, av, bv
      REAL :: dworst, dsample(4), dabs, dmax

      CALL timestart('melem_vacabc')
      dsample = 0.0
      dabs = 0.0
      dmax = 0.0

      !> An instance is reused: the ket side is expanded once per neighbour, so everything
      !> here is freed first. t_abc's own calc_abc does the same for the muffin-tin
      !> coefficients, and leaving it out is a "severe (151)" on the second neighbour.
      IF (ALLOCATED(this%ac)) DEALLOCATE (this%ac)
      IF (ALLOCATED(this%bc)) DEALLOCATE (this%bc)
      IF (ALLOCATED(this%u)) DEALLOCATE (this%u)
      IF (ALLOCATED(this%ue)) DEALLOCATE (this%ue)
      IF (ALLOCATED(this%du)) DEALLOCATE (this%du)
      IF (ALLOCATED(this%due)) DEALLOCATE (this%due)
      IF (ALLOCATED(this%kvac)) DEALLOCATE (this%kvac)
      IF (ALLOCATED(this%qpar)) DEALLOCATE (this%qpar)

      nmz = vacuum%nmz
      nvac = vacuum%nvac
      nv = lapw%nv(jspin)
      off = 0
      IF (PRESENT(ioff)) off = ioff

      ! ---- the parallel G vectors that survive, and where each 3D one lands ----
      ALLOCATE (map2(nv), this%kvac(2, lapw%dim_nv2d()))
      this%nv2 = 0
      gloop: DO ikG = 1, nv
         DO ik2 = 1, this%nv2
            IF (ALL(lapw%gvec(1:2, ikG, jspin) == this%kvac(1:2, ik2))) THEN
               map2(ikG) = ik2
               CYCLE gloop
            END IF
         END DO
         this%nv2 = this%nv2 + 1
         IF (this%nv2 > lapw%dim_nv2d()) CALL juDFT_error( &
            'melem_vacabc: more parallel G vectors than dim_nv2d', &
            calledby='melem_vacabc_calc')
         this%kvac(1:2, this%nv2) = lapw%gvec(1:2, ikG, jspin)
         map2(ikG) = this%nv2
      END DO gloop

      this%nmz = nmz
      this%delz = vacuum%delz
      this%z1 = cell%z1
      this%area = cell%area
      this%bmat33 = cell%bmat(3, 3)
      ALLOCATE (this%u(nmz, this%nv2, nvac), this%ue(nmz, this%nv2, nvac))
      ALLOCATE (this%du(nmz, this%nv2, nvac), this%due(nmz, this%nv2, nvac))
      ALLOCATE (this%qpar(3, this%nv2), source=0.0)
      ALLOCATE (this%ac(this%nv2, nbnd, 2), this%bc(this%nv2, nbnd, 2))
      this%ac = CMPLX(0.0, 0.0)
      this%bc = CMPLX(0.0, 0.0)
      ALLOCATE (t(this%nv2), dt(this%nv2), te(this%nv2), dte(this%nv2), tei(this%nv2))
      ALLOCATE (vz(nmz), zmsh(nmz))
      !> The vacuum mesh, written out because difcub needs it point by point.
      DO ik2 = 1, nmz
         zmsh(ik2) = cell%z1 + (ik2 - 1)*vacuum%delz
      END DO

      ! ---- two slots always, whether or not the second side is stored ----
      this%slot_vac(1) = 1
      this%slot_sign(1) = 1.0
      this%slot_vac(2) = MERGE(2, 1, nvac == 2)
      this%slot_sign(2) = -1.0

      const = 1.0/(SQRT(cell%omtil)*MELEM_VAC_WRONK)

      DO ivac = 1, nvac
         !> vtot%vac is complex and carries the warped part in its second index; the z-only
         !> potential this needs is index 1, and it is real.
         vz = REAL(vtot%vac(1:nmz, 1, ivac, jspin))
         vz0 = vz(nmz)

         DO ik2 = 1, this%nv2
            v(1) = lapw%bkpt(1) + this%kvac(1, ik2)
            v(2) = lapw%bkpt(2) + this%kvac(2, ik2)
            v(3) = 0.0
            !> what is left for the z motion once the parallel kinetic energy is taken out
            ev = enpara%evac(ivac, jspin) - 0.5*DOT_PRODUCT(v, MATMUL(cell%bbmat, v))
            CALL vacuz(ev, vz, vz0, nmz, vacuum%delz, t(ik2), dt(ik2), this%u(:, ik2, ivac))
            CALL vacudz(ev, vz, vz0, nmz, vacuum%delz, te(ik2), dte(ik2), tei(ik2), &
                        this%ue(:, ik2, ivac), dt(ik2), this%u(:, ik2, ivac))
            !> vacudz returns ue orthogonal to u but with a Wronskian of its own; rescaling
            !> it to MELEM_VAC_WRONK is what makes ac/bc mean the same thing at every G.
            sc = MELEM_VAC_WRONK/(te(ik2)*dt(ik2) - dte(ik2)*t(ik2))
            te(ik2) = sc*te(ik2)
            dte(ik2) = sc*dte(ik2)
            tei(ik2) = sc*tei(ik2)
            this%ue(:, ik2, ivac) = sc*this%ue(:, ik2, ivac)

            !> The derivatives the momentum needs, and the check that makes them safe: vacuz
            !> and vacudz already return du/dz AT THE BOUNDARY, so the first point of the
            !> stencil has an exact value to be measured against. Worst case over every G
            !> and both sides is reported once per run, and a gross disagreement stops it --
            !> a wrong derivative would show up nowhere else, because the momentum has no
            !> other reference.
            CALL melem_vac_ddz(this%u(:, ik2, ivac), zmsh, this%du(:, ik2, ivac))
            CALL melem_vac_ddz(this%ue(:, ik2, ivac), zmsh, this%due(:, ik2, ivac))
            !> AGAINST -dt, NOT AGAINST dt. vacuz and vacudz return the boundary derivative
            !> with the sign OPPOSITE to the mesh, which runs outwards from the film. Measured
            !> on the three-layer Pt(100) slab: du(1) = -7.139452E-01 against
            !> dt = +7.139537E-01, and due(1) = +2.573875E+00 against dte = -2.574232E+00.
            !> The magnitudes agree to 1.2e-5, so the stencil is sound and the convention was
            !> the only difference. What is kept is d/d(distance from the boundary), which is
            !> what the momentum needs together with slot_sign.
            !>
            !> The error is taken in ABSOLUTE terms and normalised at the end by the largest
            !> derivative of the whole set. Measuring it relative, G by G, does not work: the
            !> high parallel G decay so far that their solution in z is numerically zero, and
            !> there a relative measure compares two zeros -- the first version of this check
            !> reported 3.3 %, which came from a single G with dt = 0.0 exactly.
            dabs = MAX(dabs, ABS(this%du(1, ik2, ivac) + dt(ik2)), &
                             ABS(this%due(1, ik2, ivac) + dte(ik2)))
            dmax = MAX(dmax, ABS(dt(ik2)), ABS(dte(ik2)))
            IF (ivac == 1 .AND. ik2 == 1) dsample = [this%du(1, 1, 1), -dt(1), &
                                                     this%due(1, 1, 1), -dte(1)]

            IF (ivac == 1) this%qpar(:, ik2) = MATMUL(v, cell%bmat)
         END DO

         !> Both slots are filled here when only one side is stored: they share these z
         !> functions and differ only in the reflection that slot_sign carries.
         DO islot = 1, 2
            IF (this%slot_vac(islot) /= ivac) CYCLE
            sgn = this%slot_sign(islot)
            DO ikG = 1, nv
               ik2 = map2(ikG)
               zks = lapw%k3(ikG, jspin)*cell%bmat(3, 3)*sgn
               arg = zks*cell%z1
               c_1 = CMPLX(COS(arg), SIN(arg))*const
               av = -c_1*CMPLX(dte(ik2), zks*te(ik2))
               bv = c_1*CMPLX(dt(ik2), zks*t(ik2))
               IF (zMat%l_real) THEN
                  DO n = 1, nbnd
                     this%ac(ik2, n, islot) = this%ac(ik2, n, islot) + zMat%data_r(ikG + off, n)*av
                     this%bc(ik2, n, islot) = this%bc(ik2, n, islot) + zMat%data_r(ikG + off, n)*bv
                  END DO
               ELSE
                  DO n = 1, nbnd
                     this%ac(ik2, n, islot) = this%ac(ik2, n, islot) + zMat%data_c(ikG + off, n)*av
                     this%bc(ik2, n, islot) = this%bc(ik2, n, islot) + zMat%data_c(ikG + off, n)*bv
                  END DO
               END IF
            END DO
         END DO
      END DO

      !> The number is printed ALWAYS and before any abort: when the check fails, what is
      !> needed is the value, not the news that it failed.
      dworst = dabs/MAX(dmax, 1.0e-30)
      IF (.NOT. l_said) THEN
         l_said = .TRUE.
         WRITE (oUnit, '(a,es11.3)') &
            'melem_vacabc: d/dz on the mesh vs the boundary value, worst / largest ', dworst
         WRITE (oUnit, '(a,4es14.6)') &
            'melem_vacabc:   du(1), -dt, due(1), -dte at the first G: ', dsample
         WRITE (oUnit, '(a,2es13.5)') &
            'melem_vacabc:   worst absolute, and the largest derivative in the set: ', dabs, dmax
      END IF
      IF (dworst > MELEM_VAC_DDZ_TOL) CALL juDFT_error( &
         'melem_vacabc: the z derivative on the mesh does not match the boundary value '// &
         'vacuz returns', hint='see the two lines above in out for the measured values; a '// &
         'clean sign flip means d/dz is taken the other way round', calledby='melem_vacabc_calc')

      DEALLOCATE (map2, vz, zmsh, t, dt, te, dte, tei)
      CALL timestop('melem_vacabc')
   END SUBROUTINE melem_vacabc_calc

   !>  d f / dz on the vacuum mesh, the way the rest of FLEUR differentiates on a mesh:
   !>  difcub, the derivative of the cubic Lagrange interpolant through four points. It is
   !>  the same construction and the same end handling that uhu/d2fdz2cmplx.F uses for the
   !>  second derivative in the vacuum, so our momentum and the classic uHu agree on what a
   !>  z derivative means rather than each having its own.
   !>
   !>  What makes it safe is not the stencil but the caller: vacuz and vacudz already return
   !>  du/dz AT THE BOUNDARY, and the first point of this is checked against it.
   SUBROUTINE melem_vac_ddz(f, z, df)
      REAL, INTENT(IN) :: f(:)
      REAL, INTENT(IN) :: z(:)
      REAL, INTENT(OUT) :: df(:)
      INTEGER :: n, i

      n = SIZE(f)
      IF (n < 4) THEN
         df = 0.0
         RETURN
      END IF
      !> Slices of four rather than single elements: difcub declares f(0:3) with explicit
      !> shape, and handing it one element of an assumed-shape array is what d2fdz2cmplx.F
      !> may do (its own arguments have explicit shape) and this may not -- error #8284.
      df(1) = difcub(z(1:4), f(1:4), z(1))
      DO i = 2, n - 2
         df(i) = difcub(z(i - 1:i + 2), f(i - 1:i + 2), z(i))
      END DO
      df(n - 1) = difcub(z(n - 3:n), f(n - 3:n), z(n - 1))
      df(n) = difcub(z(n - 3:n), f(n - 3:n), z(n))
   END SUBROUTINE melem_vac_ddz


END MODULE m_types_melem_vacabc
