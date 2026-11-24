!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_hs_int
      use m_juDFT
      USE m_gf_types
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: gf_hs_int
      CONTAINS

      SUBROUTINE gf_hs_int(lapw, stars, napw, stepf, stepmt, jspin, fmpi,&
     &                     bkpt, bbmat, vpw, hmat, smat)
!*********************************************************************
!     Interstitial part of the Hamiltonian and overlap matrix of one
!     layer. Replaces the old hsintstep working on packed storage.
!
!     Differences to the standard FLEUR hs_int:
!     - the step function is the GF one (zero in the auxiliary
!       volumes) and is given on the FULL G-difference box
!       (-mx1:mx1, -mx2:mx2, -2*napw:2*napw) - kinetic energy and
!       overlap must use it also for differences outside the star
!       list, which the standard hs_int cannot do.
!     - the potential vpw is the GF potential (already warped and
!       zeroed in the auxiliary volumes).
!     - kinetic energy: Dirac form 1/2 (k+G_i)(k+G_j) at the embedding
!       planes, symmetrized Laplace form at the MT spheres as for the
!       LAPW MT part (hsmt):
!         T = 1/2 (k+G_i)(k+G_j) Theta_I - 1/4 |G_i-G_j|^2 Theta_MT
!
!     The matrix elements are stored in the standard modern
!     convention: data_c(jrow, icol) = <G_row|M|G_col> for row<=col
!     (upper triangle), matching hsmt which adds the MT part.
!     Noco: hmat/smat are the 2x2 spin blocks as in FLEUR's hs_int;
!     the off-diagonal blocks get V21=vpw(:,3) only.
!
!     D. Wortmann, Juelich, 2001
!*********************************************************************
      IMPLICIT NONE
!     ..
!     .. Arguments ..
      TYPE(t_lapw),INTENT(IN)       :: lapw
      TYPE(t_stars),INTENT(IN)      :: stars
      INTEGER, INTENT (IN)          :: napw
      COMPLEX, INTENT (IN)          :: stepf(-stars%mx1:,-stars%mx2:,-2 &
     &     *napw:)
      COMPLEX, INTENT (IN)          :: stepmt(-stars%mx1:,-stars%mx2:,-2&
     &     *napw:)
      INTEGER, INTENT (IN)          :: jspin
      TYPE(t_mpi), INTENT(IN)       :: fmpi
      REAL,    INTENT (IN)          :: bkpt(3),bbmat(3,3)
      COMPLEX, INTENT (IN)          :: vpw(:,:)
      CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:),smat(:,:)
!     ..
!     .. Locals ..
      INTEGER :: isp,ispPr

      IF (SIZE(hmat,1)==1) THEN
         CALL priv_diag(jspin,hmat(1,1),smat(1,1))
         RETURN
      ENDIF
      DO isp = 1,2
         CALL priv_diag(isp,hmat(isp,isp),smat(isp,isp))
      ENDDO
      !spin off-diagonal blocks (row spin ispPr, column spin isp)
      DO isp = 1,2
         ispPr = 3-isp
         CALL priv_offdiag(ispPr,isp,hmat(ispPr,isp))
      ENDDO

      CONTAINS

      SUBROUTINE priv_diag(isp,hm,sm)
      INTEGER,INTENT(IN)         :: isp
      CLASS(t_mat),INTENT(INOUT) :: hm,sm

      COMPLEX th,ts
      REAL bi(3),bj(3),r2,dg2
      INTEGER i,j,i1,i2,i3,in

      !loop over the columns (k+G_i)
      DO i = fmpi%n_rank+1, lapw%nv(isp), fmpi%n_size
         bi = bkpt + lapw%gvec(:,i,isp)
         !loop over the rows (k+G_j), j<=i
         DO j = 1, i
            !difference G_row - G_col as in the standard convention
            i1 = lapw%gvec(1,j,isp) - lapw%gvec(1,i,isp)
            i2 = lapw%gvec(2,j,isp) - lapw%gvec(2,i,isp)
            i3 = lapw%gvec(3,j,isp) - lapw%gvec(3,i,isp)
            bj = bkpt + lapw%gvec(:,j,isp)
            r2 = dot_product(bj,matmul(bi,bbmat))
            dg2 = dot_product(bj-bi,matmul(bj-bi,bbmat))
            !kinetic energy and overlap on the full box
            ts = stepf(i1,i2,i3)
            th = 0.5*r2*ts - 0.25*dg2*stepmt(i1,i2,i3)
            !add the potential if there is a corresponding star
            IF (ABS(i1)<=stars%mx1.AND.ABS(i2)<=stars%mx2.AND.           &
     &          ABS(i3)<=stars%mx3) THEN
               in = stars%ig(i1,i2,i3)
               IF (in>0) th = th + vpw(in,isp)
            ENDIF
            hm%data_c(j,i) = th
            sm%data_c(j,i) = ts
         ENDDO
      ENDDO
      END SUBROUTINE priv_diag

      SUBROUTINE priv_offdiag(ispPr,isp,hm)
      !<ispPr,G_j|V|isp,G_i> for j<=i: V21(G_j-G_i) for (2,1),
      !conjg(V21(G_i-G_j)) for (1,2)
      INTEGER,INTENT(IN)         :: ispPr,isp
      CLASS(t_mat),INTENT(INOUT) :: hm

      INTEGER i,j,g(3),in,fact

      fact = MERGE(1,-1,ispPr==2)
      DO i = fmpi%n_rank+1, lapw%nv(isp), fmpi%n_size
         DO j = 1, MIN(i,lapw%nv(ispPr))
            g = fact*(lapw%gvec(:,j,ispPr) - lapw%gvec(:,i,isp))
            IF (ABS(g(1))>stars%mx1.OR.ABS(g(2))>stars%mx2.OR.           &
     &          ABS(g(3))>stars%mx3) CYCLE
            in = stars%ig(g(1),g(2),g(3))
            IF (in==0) CYCLE
            IF (ispPr==2) THEN
               hm%data_c(j,i) = vpw(in,3)
            ELSE
               hm%data_c(j,i) = CONJG(vpw(in,3))
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE priv_offdiag

      END SUBROUTINE gf_hs_int
      END
