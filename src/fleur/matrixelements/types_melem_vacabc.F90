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
!>  momentum both need them, which is why they are built once here rather than inside each
!>  consumer: the classic tree builds them twice, inline in wann_mmkb_vac and again in
!>  wann_2dvacabcof for the plots.
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
   USE m_types
   USE m_vacuz
   USE m_vacudz
   IMPLICIT NONE
   PRIVATE

   !> The Wronskian the two z solutions are rescaled to. It is what fixes the meaning of
   !> ac/bc, so it belongs to the type and not to a caller.
   REAL, PARAMETER, PUBLIC :: MELEM_VAC_WRONK = 2.0

   TYPE, PUBLIC :: t_melem_vacabc
      COMPLEX, ALLOCATABLE :: ac(:, :, :)   !< (nv2, nbnd, 2) resolved by slot
      COMPLEX, ALLOCATABLE :: bc(:, :, :)   !< (nv2, nbnd, 2)
      REAL, ALLOCATABLE :: u(:, :, :)       !< (nmz, nv2, nvac) the decaying solution
      REAL, ALLOCATABLE :: ue(:, :, :)      !< (nmz, nv2, nvac) its energy derivative
      INTEGER, ALLOCATABLE :: kvac(:, :)    !< (2, nv2) the distinct parallel G
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
      REAL, ALLOCATABLE :: vz(:), t(:), dt(:), te(:), dte(:), tei(:)
      COMPLEX :: c_1, av, bv

      CALL timestart('melem_vacabc')

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
      ALLOCATE (this%ac(this%nv2, nbnd, 2), this%bc(this%nv2, nbnd, 2))
      this%ac = CMPLX(0.0, 0.0)
      this%bc = CMPLX(0.0, 0.0)
      ALLOCATE (t(this%nv2), dt(this%nv2), te(this%nv2), dte(this%nv2), tei(this%nv2))
      ALLOCATE (vz(nmz))

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

      DEALLOCATE (map2, vz, t, dt, te, dte, tei)
      CALL timestop('melem_vacabc')
   END SUBROUTINE melem_vacabc_calc

END MODULE m_types_melem_vacabc
