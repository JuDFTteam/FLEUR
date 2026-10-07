!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_vacbasis
   !! The LAPW basis as seen from the vacuum region of a film.
   !!
   !! Bundles the three things every vacuum routine needs for one vacuum `ivac`:
   !!
   !!  * the 3D -> 2D LAPW mapping (`nv2`, `kvac`, `map2`),
   !!  * the vacuum radial functions `u`, `ue` and their boundary values
   !!    `uz`, `duz`, `udz`, `dudz`, `ddnv` (Wronskian-normalised),
   !!  * the A/B expansion coefficients `ac`, `bc` of the eigenvectors.
   !!
   !! Field names follow those of `vacfun`, which solves the same radial problem
   !! for the Hamiltonian setup. The work is done by the routines of `m_vac_abcof`,
   !! which the hybrid film code shares.
   !!
   !! In a DFPT calculation a second instance holds the same quantities for the
   !! shifted basis at k+q, so no separate set of "q" arrays is needed.
   USE m_juDFT
   USE m_vac_abcof
   USE m_types_cell
   USE m_types_lapw
   USE m_types_mat
   USE m_types_vacuum
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_vacbasis

   TYPE t_vacbasis
      REAL                 :: wronk = 2.0    !! Wronskian of the two radial solutions (fixed in vac_uz)
      INTEGER, ALLOCATABLE :: nv2(:)         !! (jspins) number of distinct 2D G-vectors
      INTEGER, ALLOCATABLE :: kvac(:,:,:)    !! (2, nv2d, jspins) the distinct 2D G-vectors
      INTEGER, ALLOCATABLE :: map2(:,:)      !! (nvd, jspins) 3D LAPW index -> 2D index
      REAL,    ALLOCATABLE :: u(:,:,:)       !! (nmzd, nv2d, jspins) regular solution
      REAL,    ALLOCATABLE :: ue(:,:,:)      !! (nmzd, nv2d, jspins) energy derivative
      REAL,    ALLOCATABLE :: uz(:,:)        !! (nv2d, jspins) u at the vacuum boundary
      REAL,    ALLOCATABLE :: duz(:,:)       !! (nv2d, jspins) du/dz at the boundary
      REAL,    ALLOCATABLE :: udz(:,:)       !! (nv2d, jspins) ue at the boundary
      REAL,    ALLOCATABLE :: dudz(:,:)      !! (nv2d, jspins) due/dz at the boundary
      REAL,    ALLOCATABLE :: ddnv(:,:)      !! (nv2d, jspins) norm of ue
      COMPLEX, ALLOCATABLE :: ac(:,:,:)      !! (nv2d, neig, jspins) A coefficients
      COMPLEX, ALLOCATABLE :: bc(:,:,:)      !! (nv2d, neig, jspins) B coefficients
   CONTAINS
      PROCEDURE :: init         => vacbasis_init
      PROCEDURE :: calc_radfun  => vacbasis_calc_radfun
      PROCEDURE :: calc_abcoeff => vacbasis_calc_abcoeff
   END TYPE t_vacbasis

CONTAINS

   SUBROUTINE vacbasis_init(this, lapw, jspins, nmzd, neig)
      !! Allocate the arrays and build the 3D -> 2D LAPW mapping for `lapw`.
      IMPLICIT NONE

      CLASS(t_vacbasis), INTENT(INOUT) :: this
      TYPE(t_lapw),      INTENT(IN)    :: lapw
      INTEGER,           INTENT(IN)    :: jspins, nmzd, neig

      INTEGER :: nv2d, nvd, jspin

      nv2d = lapw%dim_nv2d()
      nvd  = lapw%dim_nvd()

      ALLOCATE(this%nv2(jspins), this%kvac(2,nv2d,jspins), this%map2(nvd,jspins))
      ALLOCATE(this%u(nmzd,nv2d,jspins), this%ue(nmzd,nv2d,jspins))
      ALLOCATE(this%uz(nv2d,jspins), this%duz(nv2d,jspins), this%udz(nv2d,jspins), &
               this%dudz(nv2d,jspins), this%ddnv(nv2d,jspins))
      ALLOCATE(this%ac(nv2d,neig,jspins), this%bc(nv2d,neig,jspins))

      this%kvac = 0
      this%map2 = 0
      this%ac   = CMPLX(0.0,0.0)
      this%bc   = CMPLX(0.0,0.0)

      DO jspin = 1, jspins
         CALL vac_map2(lapw, jspin, nv2d, this%kvac(1,:,jspin), this%kvac(2,:,jspin), &
                       this%map2(:,jspin), this%nv2(jspin))
      END DO

   END SUBROUTINE vacbasis_init

   SUBROUTINE vacbasis_calc_radfun(this, vacuum, cell, bkpt, gshift, evac, vz, jsp_start, jsp_end)
      !! Solve the 1D vacuum Schroedinger equation for every 2D G-vector, with the
      !! solutions scaled to satisfy the Wronskian relation.
      !!
      !! `evac` and `vz` are already restricted to the vacuum being treated;
      !! `gshift` carries the spin-spiral q/2 (non-collinear) or the phonon q
      !! (DFPT) added to `bkpt`.
      IMPLICIT NONE

      CLASS(t_vacbasis), INTENT(INOUT) :: this
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_cell),      INTENT(IN)    :: cell
      REAL,              INTENT(IN)    :: bkpt(3)
      REAL,              INTENT(IN)    :: gshift(:,:) !! (2, jspins)
      REAL,              INTENT(IN)    :: evac(:)     !! (jspins)
      REAL,              INTENT(IN)    :: vz(:,:)     !! (nmzd, jspins)
      INTEGER,           INTENT(IN)    :: jsp_start, jsp_end

      INTEGER :: jspin

      DO jspin = jsp_start, jsp_end
         CALL vac_uz(vacuum, cell, evac(jspin), vz(:,jspin), bkpt(1:2) + gshift(1:2,jspin), &
                     this%kvac(1,:,jspin), this%kvac(2,:,jspin), this%nv2(jspin), &
                     this%u(:,:,jspin), this%ue(:,:,jspin), this%uz(:,jspin), this%duz(:,jspin), &
                     this%udz(:,jspin), this%dudz(:,jspin), this%ddnv(:,jspin))
      END DO

   END SUBROUTINE vacbasis_calc_radfun

   SUBROUTINE vacbasis_calc_abcoeff(this, lapw, cell, zMat, ne, jspin, ivac, row_offset, prefac)
      !! Accumulate the A/B coefficients of the `ne` lowest eigenvectors in vacuum `ivac`.
      !!
      !! `row_offset` skips to the spin-down half of the eigenvector in the
      !! non-collinear case; `prefac` is 2 for the DFPT response coefficients.
      IMPLICIT NONE

      CLASS(t_vacbasis), INTENT(INOUT) :: this
      TYPE(t_lapw),      INTENT(IN)    :: lapw
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_mat),       INTENT(IN)    :: zMat
      INTEGER,           INTENT(IN)    :: ne, jspin, ivac, row_offset
      REAL,              INTENT(IN)    :: prefac

      CALL vac_abcof(cell, lapw, jspin, ivac, SIZE(this%ac,1), ne, row_offset, prefac, &
                     this%map2(:,jspin), this%uz(:,jspin), this%duz(:,jspin), &
                     this%udz(:,jspin), this%dudz(:,jspin), zMat, &
                     this%ac(:,:,jspin), this%bc(:,:,jspin))

   END SUBROUTINE vacbasis_calc_abcoeff

END MODULE m_types_vacbasis
