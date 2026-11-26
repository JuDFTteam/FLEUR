!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_anderson_sf
   !! Saddle-free Anderson mixing (imix="AndersonSF").
   !!
   !! Anderson mixing (like every quasi-Newton mixer) converges to any self-consistent density, also to
   !! a saddle point of the free energy. There the Jacobian of G = rho_in - rho_out has a negative
   !! eigenvalue and the Newton step moves towards the saddle. This routine
   !!  1. computes the Anderson step from its least-squares coefficients,
   !!  2. projects the Jacobian onto the mixing history, (dX^T w dG) v = lambda (dX^T w dX) v, and
   !!     accepts a mode as unstable only if lambda is clearly negative, nearly real and consistent
   !!     with the history, the residual is above the noise level of the run, and the direction was not
   !!     rejected before,
   !!  3. escapes along an unstable direction u: the Anderson step is restricted to the complement of u,
   !!     and along u a push is added in the downhill direction (residual along u positive). The push
   !!     doubles as long as the residual along u does not fall, i.e. the minimum along u is not near.
   !!     The first escape step is the push alone, so that the response of the residual along u is not
   !!     mixed with the response to the Anderson step in the complement of u,
   !!  4. ends the escape when the residual along u has fallen below 20% of its maximum (minimum along u
   !!     reached; the history is restarted), or rejects u as a false detection when the residual along
   !!     u turns against the push without ever growing,
   !!  5. limits the step outside of the escape direction to max_step_factor times alpha*|F|.
   !! Without an unstable mode and below the step limit the step is the Anderson step.
   !! Notation: dX and dF are the differences of successive input densities and residuals
   !! F = rho_out - rho_in in the history, dG = -dF, w is the mixing metric.
   USE m_juDFT
   USE m_types_mixvector
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: anderson_sf

   ! Acceptance of an unstable mode
   REAL,    PARAMETER :: neg_eig_tolerance = 0.01      !! Re(lambda) < -neg_eig_tolerance*max|lambda|
   REAL,    PARAMETER :: imag_part_tolerance = 0.5     !! |Im(lambda)| <= imag_part_tolerance*|Re(lambda)|
   REAL,    PARAMETER :: max_mode_inconsistency = 1.0  !! |dG v - lambda dX v| <= this*|lambda dX v|
   REAL,    PARAMETER :: noise_level = 1.0e-5          !! no detection for |F| < noise_level*max|F| of the run
   REAL,    PARAMETER :: rejected_overlap = 0.8        !! overlap above which a direction counts as rejected
   ! Step control
   REAL,    PARAMETER :: max_step_factor = 16.0        !! step limit in units of alpha*|F|
   INTEGER, PARAMETER :: max_push_doublings = 4
   INTEGER, PARAMETER :: max_escape_steps = 50
   INTEGER, PARAMETER :: max_rejected_dirs = 8

   ! Escape state, kept between iterations
   LOGICAL, SAVE           :: l_escape = .FALSE.
   TYPE(t_mixvector), SAVE :: escape_dir             !! unit escape direction u
   TYPE(t_mixvector), SAVE :: w_escape_dir           !! w*u
   REAL, SAVE              :: res_along_dir_start    !! residual along u when the escape started
   REAL, SAVE              :: res_along_dir_max      !! its maximum during the escape
   REAL, SAVE              :: res_along_dir_prev     !! its value in the previous escape step
   REAL, SAVE              :: push_start             !! initial push length
   INTEGER, SAVE           :: n_escape_steps, n_push_doublings
   REAL, SAVE              :: max_res_norm = 0.0     !! largest |F| of the run, reference for the noise level
   INTEGER, SAVE           :: n_rejected = 0
   TYPE(t_mixvector), SAVE :: w_rejected_dirs(max_rejected_dirs)  !! w*u of rejected directions

CONTAINS

   SUBROUTINE anderson_sf(alpha, res, rho_in, l_dfpt, new_mode_eig, step_rel, l_restart, escape_state, &
                          push_rel, res_along_dir_rel)
      REAL, INTENT(IN)                 :: alpha
      TYPE(t_mixvector), INTENT(IN)    :: res(:)       !! residuals F = rho_out - rho_in of the history, latest last
      TYPE(t_mixvector), INTENT(INOUT) :: rho_in(:)    !! input densities of the history; the last one is updated
      LOGICAL, INTENT(IN)              :: l_dfpt
      REAL, INTENT(OUT)                :: new_mode_eig !! eigenvalue of a newly accepted unstable mode, else 0
      REAL, INTENT(OUT)                :: step_rel     !! |step|/|alpha*F| before the step limit
      LOGICAL, INTENT(OUT)             :: l_restart    !! the history has to be restarted
      INTEGER, INTENT(OUT)             :: escape_state !! 0 none, 1 escaping, 2 completed, 3 rejected, 4 aborted
      REAL, INTENT(OUT)                :: push_rel     !! push length along u in units of alpha*|F|
      REAL, INTENT(OUT)                :: res_along_dir_rel !! residual along u relative to its maximum

      TYPE(t_mixvector), ALLOCATABLE :: dx(:), df(:), w_dx(:), w_df(:)
      TYPE(t_mixvector) :: step, w_res, mode_dir
      REAL, ALLOCATABLE :: dx_dx(:,:), dx_dg(:,:), df_df(:,:), df_res(:), anderson_coef(:)
      REAL    :: res_norm, res_along_dir, step_along_dir, mode_coef_abs, mode_dir_norm
      REAL    :: step_norm, step_limit, push_length
      INTEGER :: n_hist, n_diff, i, j

      n_hist = SIZE(res)
      new_mode_eig = 0.0
      step_rel = 1.0
      l_restart = .FALSE.
      escape_state = MERGE(1, 0, l_escape)
      push_rel = 0.0
      res_along_dir_rel = 0.0

      w_res = res(n_hist)%apply_metric(l_dfpt)
      res_norm = SQRT(ABS(res(n_hist).dot.w_res))
      max_res_norm = MAX(max_res_norm, res_norm)
      IF (n_hist < 2) THEN ! no history: linear mixing
         rho_in(n_hist) = rho_in(n_hist) + alpha*res(n_hist)
         RETURN
      END IF

      ! End of an escape: false detection, minimum along u reached, or too long
      IF (l_escape) THEN
         res_along_dir = w_escape_dir.dot.res(n_hist)
         res_along_dir_max = MAX(res_along_dir_max, res_along_dir)
         res_along_dir_rel = res_along_dir/res_along_dir_max
         IF (n_escape_steps >= 2 .AND. res_along_dir <= 0.0 .AND. res_along_dir_max <= 1.05*res_along_dir_start) THEN
            n_rejected = MIN(n_rejected + 1, max_rejected_dirs)
            w_rejected_dirs(n_rejected) = w_escape_dir
            l_escape = .FALSE.
            escape_state = 3
            l_restart = .TRUE.
         ELSE IF (n_escape_steps >= 2 .AND. res_along_dir < 0.2*res_along_dir_max) THEN
            l_escape = .FALSE.
            escape_state = 2
            l_restart = .TRUE.
         ELSE IF (n_escape_steps >= max_escape_steps) THEN
            l_escape = .FALSE.
            escape_state = 4
         END IF
         IF (l_restart) THEN
            rho_in(n_hist) = rho_in(n_hist) + alpha*res(n_hist)
            RETURN
         END IF
      END IF

      ! Differences of the history and their scalar products in the metric w
      n_diff = n_hist - 1
      ALLOCATE(dx(n_diff), df(n_diff), w_dx(n_diff), w_df(n_diff))
      DO i = 1, n_diff
         dx(i) = rho_in(i+1) - rho_in(i)
         df(i) = res(i+1) - res(i)
         w_dx(i) = dx(i)%apply_metric(l_dfpt)
         w_df(i) = df(i)%apply_metric(l_dfpt)
      END DO
      ALLOCATE(dx_dx(n_diff,n_diff), dx_dg(n_diff,n_diff), df_df(n_diff,n_diff), df_res(n_diff), anderson_coef(n_diff))
      DO j = 1, n_diff
         DO i = 1, n_diff
            dx_dx(i,j) = w_dx(i).dot.dx(j)
            dx_dg(i,j) = -(w_dx(i).dot.df(j))
            df_df(i,j) = w_df(i).dot.df(j)
         END DO
         df_res(j) = w_df(j).dot.res(n_hist)
      END DO

      ! Anderson step x_new = x + alpha*F - sum_i g_i (dx_i + alpha*df_i), (dF^T w dF) g = dF^T w F
      step = alpha*res(n_hist)
      IF (solve_linear(df_df, df_res, anderson_coef)) THEN
         DO i = 1, n_diff
            step = step - anderson_coef(i)*(dx(i) + alpha*df(i))
         END DO

         ! Start of an escape along a newly detected unstable mode
         IF (.NOT.l_escape .AND. res_norm >= noise_level*max_res_norm) THEN
            IF (find_unstable_mode(dx, df, dx_dx, dx_dg, anderson_coef, l_dfpt, mode_dir, new_mode_eig, mode_coef_abs)) THEN
               mode_dir_norm = SQRT(ABS(mode_dir.dot.mode_dir%apply_metric(l_dfpt)))
               escape_dir = (1.0/mode_dir_norm)*mode_dir
               w_escape_dir = escape_dir%apply_metric(l_dfpt)
               res_along_dir = w_escape_dir.dot.res(n_hist)
               IF (res_along_dir < 0.0) THEN ! downhill direction
                  escape_dir = (-1.0)*escape_dir
                  w_escape_dir = (-1.0)*w_escape_dir
                  res_along_dir = -res_along_dir
               END IF
               l_escape = .TRUE.
               escape_state = 1
               res_along_dir_start = res_along_dir
               res_along_dir_max = res_along_dir
               res_along_dir_prev = res_along_dir
               ! initial push: the reversed Newton component of the mode
               push_start = MAX(2.0*mode_coef_abs*mode_dir_norm, alpha*res_norm)
               n_escape_steps = 0
               n_push_doublings = 0
            END IF
         END IF
      END IF

      step_limit = max_step_factor*alpha*res_norm
      IF (l_escape) THEN
         ! Anderson in the complement of u, push along u
         step_along_dir = w_escape_dir.dot.step
         step = step - step_along_dir*escape_dir
         step_norm = SQRT(ABS(step.dot.step%apply_metric(l_dfpt)))
         IF (alpha*res_norm > 0.0) step_rel = step_norm/(alpha*res_norm)
         IF (step_norm > step_limit) step = (step_limit/step_norm)*step
         IF (n_escape_steps == 0) step = 0.0*step ! first escape step: push alone
         res_along_dir = w_escape_dir.dot.res(n_hist)
         IF (n_escape_steps > 0 .AND. res_along_dir >= 0.9*res_along_dir_prev) n_push_doublings = n_push_doublings + 1
         push_length = push_start*2.0**MIN(n_push_doublings, max_push_doublings)
         step = step + push_length*escape_dir
         res_along_dir_prev = res_along_dir
         n_escape_steps = n_escape_steps + 1
         IF (alpha*res_norm > 0.0) push_rel = push_length/(alpha*res_norm)
         res_along_dir_rel = res_along_dir/res_along_dir_max
      ELSE
         step_norm = SQRT(ABS(step.dot.step%apply_metric(l_dfpt)))
         IF (alpha*res_norm > 0.0) step_rel = step_norm/(alpha*res_norm)
         IF (step_norm > step_limit) step = (step_limit/step_norm)*step
      END IF
      rho_in(n_hist) = rho_in(n_hist) + step
   END SUBROUTINE anderson_sf

   LOGICAL FUNCTION find_unstable_mode(dx, df, dx_dx, dx_dg, anderson_coef, l_dfpt, mode_dir, mode_eig, mode_coef_abs)
      !! Most negative acceptable eigenmode of the projected Jacobian T = (dX^T w dX)^-1 dX^T w dG.
      !! mode_dir = dX Re(v); mode_coef_abs = |c| with the Anderson coefficients g = V c in the eigenbasis.
      TYPE(t_mixvector), INTENT(IN)  :: dx(:), df(:)
      REAL, INTENT(IN)               :: dx_dx(:,:), dx_dg(:,:), anderson_coef(:)
      LOGICAL, INTENT(IN)            :: l_dfpt
      TYPE(t_mixvector), INTENT(OUT) :: mode_dir
      REAL, INTENT(OUT)              :: mode_eig, mode_coef_abs

      TYPE(t_mixvector) :: mismatch
      REAL, ALLOCATABLE    :: jacobian(:,:), eig_re(:), eig_im(:), eigvec_lapack(:,:), work(:), eigvec_re(:)
      COMPLEX, ALLOCATABLE :: eigvec(:,:), eigvec_lu(:,:), mode_coef(:,:)
      INTEGER, ALLOCATABLE :: ipiv(:), order(:)
      REAL    :: dummy(1,1), dir_norm, mismatch_norm, eig_abs_max
      INTEGER :: n, i, j, jj, info, tmp, irej
      LOGICAL :: l_rejected

      find_unstable_mode = .FALSE.
      mode_eig = 0.0
      mode_coef_abs = 0.0
      n = SIZE(anderson_coef)
      ALLOCATE(jacobian(n,n), eig_re(n), eig_im(n), eigvec_lapack(n,n), work(8*n), eigvec(n,n), eigvec_lu(n,n), &
               mode_coef(n,1), ipiv(n), order(n), eigvec_re(n))
      IF (.NOT.solve_linear_multi(dx_dx, dx_dg, jacobian)) RETURN
      CALL dgeev('N', 'V', n, jacobian, n, eig_re, eig_im, dummy, 1, eigvec_lapack, n, work, SIZE(work), info)
      IF (info /= 0) RETURN
      eig_abs_max = MAXVAL(SQRT(eig_re**2 + eig_im**2))
      IF (MINVAL(eig_re) >= -neg_eig_tolerance*eig_abs_max) RETURN

      ! complex eigenvectors from the dgeev storage of conjugate pairs
      j = 1
      DO WHILE (j <= n)
         IF (eig_im(j) == 0.0) THEN
            eigvec(:,j) = CMPLX(eigvec_lapack(:,j), 0.0)
            j = j + 1
         ELSE
            eigvec(:,j) = CMPLX(eigvec_lapack(:,j), eigvec_lapack(:,j+1))
            eigvec(:,j+1) = CONJG(eigvec(:,j))
            j = j + 2
         END IF
      END DO
      ! Anderson coefficients in the eigenbasis
      eigvec_lu = eigvec
      mode_coef(:,1) = CMPLX(anderson_coef, 0.0)
      CALL zgesv(n, 1, eigvec_lu, n, ipiv, mode_coef, n, info)
      IF (info /= 0) RETURN

      ! candidates in the order of increasing Re(lambda)
      order = [(i, i=1, n)]
      DO i = 1, n - 1
         DO jj = i + 1, n
            IF (eig_re(order(jj)) < eig_re(order(i))) THEN
               tmp = order(i)
               order(i) = order(jj)
               order(jj) = tmp
            END IF
         END DO
      END DO
      DO i = 1, n
         j = order(i)
         IF (eig_re(j) >= -neg_eig_tolerance*eig_abs_max) EXIT
         IF (ABS(eig_im(j)) > imag_part_tolerance*ABS(eig_re(j))) CYCLE
         ! mode_dir = dX v, mismatch = dG v - lambda dX v
         eigvec_re = REAL(eigvec(:,j))
         mode_dir = eigvec_re(1)*dx(1)
         mismatch = (-eigvec_re(1))*df(1)
         DO jj = 2, n
            mode_dir = mode_dir + eigvec_re(jj)*dx(jj)
            mismatch = mismatch - eigvec_re(jj)*df(jj)
         END DO
         dir_norm = SQRT(ABS(mode_dir.dot.mode_dir%apply_metric(l_dfpt)))
         IF (dir_norm == 0.0) CYCLE
         mismatch = mismatch - eig_re(j)*mode_dir
         mismatch_norm = SQRT(ABS(mismatch.dot.mismatch%apply_metric(l_dfpt)))
         IF (mismatch_norm > max_mode_inconsistency*ABS(eig_re(j))*dir_norm) CYCLE
         l_rejected = .FALSE.
         DO irej = 1, n_rejected
            IF (ABS(w_rejected_dirs(irej).dot.mode_dir)/dir_norm > rejected_overlap) l_rejected = .TRUE.
         END DO
         IF (l_rejected) CYCLE
         mode_eig = eig_re(j)
         mode_coef_abs = ABS(mode_coef(j,1))
         find_unstable_mode = .TRUE.
         RETURN
      END DO
   END FUNCTION find_unstable_mode

   LOGICAL FUNCTION solve_linear(a, b, x)
      !! Solves a x = b; .FALSE. if a is singular
      REAL, INTENT(IN)  :: a(:,:), b(:)
      REAL, INTENT(OUT) :: x(:)
      REAL, ALLOCATABLE :: a_lu(:,:), rhs(:,:)
      INTEGER, ALLOCATABLE :: ipiv(:)
      INTEGER :: n, info

      n = SIZE(b)
      ALLOCATE(a_lu(n,n), rhs(n,1), ipiv(n))
      a_lu = a
      rhs(:,1) = b
      CALL dgesv(n, 1, a_lu, n, ipiv, rhs, n, info)
      solve_linear = info == 0
      x = rhs(:,1)
   END FUNCTION solve_linear

   LOGICAL FUNCTION solve_linear_multi(a, b, x)
      !! Solves a x = b for several right-hand sides; .FALSE. if a is singular
      REAL, INTENT(IN)  :: a(:,:), b(:,:)
      REAL, INTENT(OUT) :: x(:,:)
      REAL, ALLOCATABLE :: a_lu(:,:)
      INTEGER, ALLOCATABLE :: ipiv(:)
      INTEGER :: n, info

      n = SIZE(a,1)
      ALLOCATE(a_lu(n,n), ipiv(n))
      a_lu = a
      x = b
      CALL dgesv(n, SIZE(b,2), a_lu, n, ipiv, x, n, info)
      solve_linear_multi = info == 0
   END FUNCTION solve_linear_multi

END MODULE m_anderson_sf
