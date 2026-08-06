!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vac_map2
   !! Mapping of the 3D LAPW basis onto the distinct 2D (G_x,G_y) vectors used in
   !! the vacuum region of a film.
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: vac_map2

CONTAINS

   SUBROUTINE vac_map2(lapw, jspins, nv2, kvac, map2)
      !! Collect the distinct in-plane components of the LAPW G-vectors.
      !!
      !! For every spin, `kvac(:,ikG2,jspin)` lists the distinct (G_x,G_y) pairs
      !! occurring in `lapw%gvec`, `nv2(jspin)` counts them, and `map2(ikG,jspin)`
      !! maps each 3D basis function onto its 2D partner.
      USE m_types_lapw
      IMPLICIT NONE

      TYPE(t_lapw), INTENT(IN)  :: lapw
      INTEGER,      INTENT(IN)  :: jspins
      INTEGER,      INTENT(OUT) :: nv2(:)      !! (jspins)
      INTEGER,      INTENT(OUT) :: kvac(:,:,:) !! (2, lapw%dim_nv2d(), jspins)
      INTEGER,      INTENT(OUT) :: map2(:,:)   !! (lapw%dim_nvd(), jspins)

      INTEGER :: jspin, ikG, ikG2

      DO jspin = 1, jspins
         nv2(jspin) = 0
         k_loop: DO ikG = 1, lapw%nv(jspin)
            DO ikG2 = 1, nv2(jspin)
               IF (ALL(lapw%gvec(1:2,ikG,jspin) == kvac(1:2,ikG2,jspin))) THEN
                  map2(ikG,jspin) = ikG2
                  CYCLE k_loop
               END IF
            END DO
            nv2(jspin) = nv2(jspin) + 1
            IF (nv2(jspin) > lapw%dim_nv2d()) CALL juDFT_error("lapw%dim_nv2d() too small", calledby="vac_map2")
            kvac(1:2,nv2(jspin),jspin) = lapw%gvec(1:2,ikG,jspin)
            map2(ikG,jspin) = nv2(jspin)
         END DO k_loop
      END DO

   END SUBROUTINE vac_map2

END MODULE m_vac_map2
