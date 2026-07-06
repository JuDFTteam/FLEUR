!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_hs_int
      use m_juDFT
      IMPLICIT NONE
      PUBLIC :: gf_hs_int
      CONTAINS

      SUBROUTINE gf_hs_int(lapw, stars, napw, stepf, jspin, fmpi,        &
     &                     bkpt, bbmat, vpw, hmat, smat)
!*********************************************************************
!     Interstitial part of the Hamiltonian and overlap matrix of one
!     layer. Replaces the old hsintstep working on packed storage.
!
!     Two differences to the standard FLEUR hs_int:
!     - the step function is the GF one (zero in the auxiliary
!       volumes) and is given on the FULL G-difference box
!       (-mx1:mx1, -mx2:mx2, -2*napw:2*napw) - kinetic energy and
!       overlap must use it also for differences outside the star
!       list, which the standard hs_int cannot do.
!     - the potential vpw is the GF potential (already warped and
!       zeroed in the auxiliary volumes).
!
!     The matrix elements are stored in the standard modern
!     convention: data_c(jrow, icol) = <G_row|M|G_col> for row<=col
!     (upper triangle), matching hsmt which adds the MT part.
!
!     D. Wortmann, Juelich, 2001
!*********************************************************************
      USE m_gf_types
      IMPLICIT NONE
!     ..
!     .. Arguments ..
      TYPE(t_lapw),INTENT(IN)       :: lapw
      TYPE(t_stars),INTENT(IN)      :: stars
      INTEGER, INTENT (IN)          :: napw
      COMPLEX, INTENT (IN)          :: stepf(-stars%mx1:,-stars%mx2:,-2 &
     &     *napw:)
      INTEGER, INTENT (IN)          :: jspin
      TYPE(t_mpi), INTENT(IN)       :: fmpi
      REAL,    INTENT (IN)          :: bkpt(3),bbmat(3,3)
      COMPLEX, INTENT (IN)          :: vpw(:,:)
      CLASS(t_mat),INTENT(INOUT)    :: hmat,smat
!     ..
!     .. Locals ..
      COMPLEX th,ts
      REAL bi(3),bj(3),r2
      INTEGER i,j,i1,i2,i3,in

      !loop over the columns (k+G_i)
      DO i = fmpi%n_rank+1, lapw%nv(jspin), fmpi%n_size
         bi = bkpt + lapw%gvec(:,i,jspin)
         !loop over the rows (k+G_j), j<=i
         DO j = 1, i
            !difference G_row - G_col as in the standard convention
            i1 = lapw%gvec(1,j,jspin) - lapw%gvec(1,i,jspin)
            i2 = lapw%gvec(2,j,jspin) - lapw%gvec(2,i,jspin)
            i3 = lapw%gvec(3,j,jspin) - lapw%gvec(3,i,jspin)
            bj = bkpt + lapw%gvec(:,j,jspin)
            r2 = dot_product(bj,matmul(bi,bbmat))
            !kinetic energy and overlap on the full box
            ts = stepf(i1,i2,i3)
            th = 0.5*r2*ts
            !add the potential if there is a corresponding star
            IF (ABS(i1)<=stars%mx1.AND.ABS(i2)<=stars%mx2.AND.           &
     &          ABS(i3)<=stars%mx3) THEN
               in = stars%ig(i1,i2,i3)
               IF (in>0) th = th + vpw(in,jspin)
            ENDIF
            hmat%data_c(j,i) = th
            smat%data_c(j,i) = ts
         ENDDO
      ENDDO

      END SUBROUTINE gf_hs_int
      END
