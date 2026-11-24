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
   !!  2. once |F| has fallen below detection_fraction of its maximum in the run, projects the Jacobian
   !!     onto the mixing history, (dX^T w dG) v = lambda (dX^T w dX) v, and accepts a mode as unstable
   !!     only if lambda is clearly negative, nearly real and consistent with the history, the residual
   !!     is above the noise level of the run, and the direction was not rejected recently. Far from
   !!     convergence the history does not describe the neighbourhood of a solution, and an escape
   !!     started there selects an arbitrary state,
   !!  3. escapes along an unstable direction u: the Anderson step is restricted to the complement of u,
   !!     and along u the saddle-free Newton step F_u/|lambda_u| is taken (F_u: residual along u,
   !!     lambda_u: curvature along u, updated by secants of the escape steps), limited to max_push_growth
   !!     times the previous push. The step leads away from the saddle while lambda_u < 0 and towards
   !!     the minimum along u once lambda_u > 0,
   !!  4. ends the escape when F_u has fallen below completion_fraction of its maximum after a confirmed
   !!     negative curvature (minimum along u reached), rejects u as a false detection when the first
   !!     escape step shows a positive curvature, and aborts it after max_escape_steps. Rejected and
   !!     aborted directions are not accepted again for rejection_lifetime mixing steps. The end of an
   !!     escape restarts the history.
   !! Without an unstable mode the step is the Anderson step.
   !! Notation: dX and dF are the differences of successive input densities and residuals
   !! F = rho_out - rho_in in the history, dG = -dF, w is the mixing metric.
   USE m_juDFT
   USE m_types_mixvector
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: anderson_sf

   ! Acceptance of an unstable mode
   REAL,    PARAMETER :: detection_fraction = 0.1      !! detection only for |F| < detection_fraction*max|F|
   REAL,    PARAMETER :: neg_eig_tolerance = 0.01      !! Re(lambda) < -neg_eig_tolerance*max|lambda|
   REAL,    PARAMETER :: imag_part_tolerance = 0.5     !! |Im(lambda)| <= imag_part_tolerance*|Re(lambda)|
   REAL,    PARAMETER :: max_mode_inconsistency = 1.0  !! |dG v - lambda dX v| <= this*|lambda dX v|
   REAL,    PARAMETER :: noise_level = 1.0e-5          !! no detection for |F| < noise_level*max|F| of the run
   REAL,    PARAMETER :: rejected_overlap = 0.8        !! overlap above which a direction counts as rejected
   INTEGER, PARAMETER :: rejection_lifetime = 10       !! number of mixing steps a direction stays rejected
   INTEGER, PARAMETER :: max_rejected_dirs = 8
   ! Escape control
   REAL,    PARAMETER :: min_eig_fraction = 0.1        !! |lambda_u| >= this*|lambda| of the detection
   REAL,    PARAMETER :: max_push_growth = 2.0         !! growth limit of the push per escape step
   REAL,    PARAMETER :: completion_fraction = 0.2     !! minimum reached when F_u < this*max(F_u)
   INTEGER, PARAMETER :: max_escape_steps = 50

   ! Escape state, kept between iterations
   LOGICAL, SAVE           :: l_escape = .FALSE.
   LOGICAL, SAVE           :: l_confirmed            !! negative curvature along u measured
   TYPE(t_mixvector), SAVE :: escape_dir             !! unit escape direction u
   TYPE(t_mixvector), SAVE :: w_escape_dir           !! w*u
   REAL, SAVE              :: eig_detected           !! lambda of the detected mode
   REAL, SAVE              :: eig_along_dir          !! current curvature estimate lambda_u
   REAL, SAVE              :: res_along_dir_max      !! maximum of F_u during the escape
   REAL, SAVE              :: res_along_dir_prev     !! F_u in the previous escape step
   REAL, SAVE              :: coord_prev             !! coordinate of the input density along u, previous step
   REAL, SAVE              :: push_prev              !! push of the previous escape step
   INTEGER, SAVE           :: n_escape_steps
   INTEGER, SAVE           :: n_mix_steps = 0        !! number of calls, clock for the rejected directions
   REAL, SAVE              :: max_res_norm = 0.0     !! largest |F| of the run, reference for the detection
   TYPE(t_mixvector), SAVE :: w_rejected_dirs(max_rejected_dirs)  !! w*u of rejected directions
   INTEGER, SAVE           :: rejection_step(max_rejected_dirs) = -rejection_lifetime - 1 !! step of the rejection

CONTAINS

   SUBROUTINE anderson_sf(alpha, res, rho_in, l_dfpt, mode_eig, step_rel, l_restart, escape_state, &
                          push_rel, res_along_dir_rel)
      REAL, INTENT(IN)                 :: alpha
      TYPE(t_mixvector), INTENT(IN)    :: res(:)       !! residuals F = rho_out - rho_in of the history, latest last
      TYPE(t_mixvector), INTENT(INOUT) :: rho_in(:)    !! input densities of the history; the last one is updated
      LOGICAL, INTENT(IN)              :: l_dfpt
      REAL, INTENT(OUT)                :: mode_eig     !! curvature along u (detected or measured), else 0
      REAL, INTENT(OUT)                :: step_rel     !! |step|/|alpha*F| outside of the escape direction
      LOGICAL, INTENT(OUT)             :: l_restart    !! the history has to be restarted
      INTEGER, INTENT(OUT)             :: escape_state !! 0 none, 1 escaping, 2 completed, 3 rejected, 4 aborted
      REAL, INTENT(OUT)                :: push_rel     !! push length along u in units of alpha*|F|
      REAL, INTENT(OUT)                :: res_along_dir_rel !! F_u relative to its maximum

      TYPE(t_mixvector), ALLOCATABLE :: dx(:), df(:), w_dx(:), w_df(:)
      TYPE(t_mixvector) :: step, w_res, mode_dir
      REAL, ALLOCATABLE :: dx_dx(:,:), dx_dg(:,:), df_df(:,:), df_res(:), anderson_coef(:)
      REAL    :: res_norm, res_along_dir, coord, d_coord, mode_dir_norm
      REAL    :: step_norm, push_length
      INTEGER :: n_hist, n_diff, i, j

      n_hist = SIZE(res)
      n_mix_steps = n_mix_steps + 1
      mode_eig = 0.0
      step_rel = 1.0
      l_restart = .FALSE.
      escape_state = MERGE(1, 0, l_escape)
      push_rel = 0.0
      res_along_dir_rel = 0.0

      w_res = res(n_hist)%apply_metric(l_dfpt)
      res_norm = SQRT(ABS(res(n_hist).dot.w_res))
      max_res_norm = MAX(max_res_norm, res_norm)

      ! Escape: curvature along u from the last escape step, end of the escape
      IF (l_escape) THEN
         res_along_dir = w_escape_dir.dot.res(n_hist)
         coord = w_escape_dir.dot.rho_in(n_hist)
         res_along_dir_max = MAX(res_along_dir_max, res_along_dir)
         res_along_dir_rel = res_along_dir/res_along_dir_max
         IF (n_escape_steps > 0) THEN
            d_coord = coord - coord_prev
            IF (d_coord /= 0.0) THEN
               eig_along_dir = -(res_along_dir - res_along_dir_prev)/d_coord
               IF (eig_along_dir < 0.0) l_confirmed = .TRUE.
            END IF
         END IF
         mode_eig = eig_along_dir
         IF (n_escape_steps > 0 .AND. .NOT.l_confirmed) THEN
            escape_state = 3 ! false detection
            CALL add_rejected_dir(w_escape_dir)
         ELSE IF (l_confirmed .AND. res_along_dir < completion_fraction*res_along_dir_max) THEN
            escape_state = 2 ! minimum along u reached
         ELSE IF (n_escape_steps >= max_escape_steps) THEN
            escape_state = 4
            CALL add_rejected_dir(w_escape_dir)
         END IF
         IF (escape_state > 1) THEN
            l_escape = .FALSE.
            l_restart = .TRUE.
            rho_in(n_hist) = rho_in(n_hist) + alpha*res(n_hist)
            RETURN
         END IF
      END IF

      step = alpha*res(n_hist)
      IF (n_hist >= 2) THEN
         ! Differences of the history and their scalar products in the metric w
         n_diff = n_hist - 1
         ALLOCATE(dx(n_diff), df(n_diff), w_dx(n_diff), w_df(n_diff))
         DO i = 1, n_diff
            dx(i) = rho_in(i+1) - rho_in(i)
            df(i) = res(i+1) - res(i)
            w_dx(i) = dx(i)%apply_metric(l_dfpt)
            w_df(i) = df(i)%apply_metric(l_dfpt)
         END DO
         ALLOCATE(dx_dx(n_diff,n_diff), dx_dg(n_diff,n_diff), df_df(n_diff,n_diff), df_res(n_diff), &
                  anderson_coef(n_diff))
         DO j = 1, n_diff
            DO i = 1, n_diff
               dx_dx(i,j) = w_dx(i).dot.dx(j)
               dx_dg(i,j) = -(w_dx(i).dot.df(j))
               df_df(i,j) = w_df(i).dot.df(j)
            END DO
            df_res(j) = w_df(j).dot.res(n_hist)
         END DO

         ! Anderson step x_new = x + alpha*F - sum_i g_i (dx_i + alpha*df_i), (dF^T w dF) g = dF^T w F
         IF (solve_linear(df_df, df_res, anderson_coef)) THEN
            DO i = 1, n_diff
               step = step - anderson_coef(i)*(dx(i) + alpha*df(i))
            END DO

            ! Start of an escape along a newly detected unstable mode
            IF (.NOT.l_escape .AND. res_norm >= noise_level*max_res_norm .AND. &
                res_norm < detection_fraction*max_res_norm) THEN
               IF (find_unstable_mode(dx, df, dx_dx, dx_dg, l_dfpt, mode_dir, mode_eig)) THEN
                  mode_dir_norm = SQRT(ABS(mode_dir.dot.mode_dir%apply_metric(l_dfpt)))
                  escape_dir = (1.0/mode_dir_norm)*mode_dir
                  w_escape_dir = escape_dir%apply_metric(l_dfpt)
                  res_along_dir = w_escape_dir.dot.res(n_hist)
                  IF (res_along_dir < 0.0) THEN ! downhill direction
                     escape_dir = (-1.0)*escape_dir
                     w_escape_dir = (-1.0)*w_escape_dir
                     res_along_dir = -res_along_dir
                  END IF
                  coord = w_escape_dir.dot.rho_in(n_hist)
                  l_escape = .TRUE.
                  l_confirmed = .FALSE.
                  escape_state = 1
                  eig_detected = mode_eig
                  eig_along_dir = mode_eig
                  res_along_dir_max = res_along_dir
                  res_along_dir_rel = 1.0
                  push_prev = 0.0
                  n_escape_steps = 0
               END IF
            END IF
         END IF
      END IF

      IF (l_escape) THEN
         ! Anderson in the complement of u, saddle-free Newton step along u
         step = step - (w_escape_dir.dot.step)*escape_dir
         step_norm = SQRT(ABS(step.dot.step%apply_metric(l_dfpt)))
         IF (alpha*res_norm > 0.0) step_rel = step_norm/(alpha*res_norm)
         push_length = res_along_dir/MAX(ABS(eig_along_dir), min_eig_fraction*ABS(eig_detected))
         IF (push_prev > 0.0) push_length = MIN(push_length, max_push_growth*push_prev)
         step = step + push_length*escape_dir
         push_prev = push_length
         res_along_dir_prev = res_along_dir
         coord_prev = coord
         n_escape_steps = n_escape_steps + 1
         IF (alpha*res_norm > 0.0) push_rel = push_length/(alpha*res_norm)
      ELSE
         step_norm = SQRT(ABS(step.dot.step%apply_metric(l_dfpt)))
         IF (alpha*res_norm > 0.0) step_rel = step_norm/(alpha*res_norm)
      END IF
      rho_in(n_hist) = rho_in(n_hist) + step
   END SUBROUTINE anderson_sf

   SUBROUTINE add_rejected_dir(w_dir)
      !! Stores w*u of a rejected direction in the slot of the oldest rejection
      TYPE(t_mixvector), INTENT(IN) :: w_dir
      INTEGER :: slot

      slot = MINLOC(rejection_step, DIM=1)
      w_rejected_dirs(slot) = w_dir
      rejection_step(slot) = n_mix_steps
   END SUBROUTINE add_rejected_dir

   LOGICAL FUNCTION find_unstable_mode(dx, df, dx_dx, dx_dg, l_dfpt, mode_dir, mode_eig)
      !! Most negative acceptable eigenmode of the projected Jacobian T = (dX^T w dX)^-1 dX^T w dG.
      !! mode_dir = dX Re(v).
      TYPE(t_mixvector), INTENT(IN)  :: dx(:), df(:)
      REAL, INTENT(IN)               :: dx_dx(:,:), dx_dg(:,:)
      LOGICAL, INTENT(IN)            :: l_dfpt
      TYPE(t_mixvector), INTENT(OUT) :: mode_dir
      REAL, INTENT(OUT)              :: mode_eig

      TYPE(t_mixvector) :: mismatch
      REAL, ALLOCATABLE :: jacobian(:,:), eig_re(:), eig_im(:), eigvec(:,:), work(:)
      INTEGER, ALLOCATABLE :: order(:)
      REAL    :: dummy(1,1), dir_norm, mismatch_norm, eig_abs_max
      INTEGER :: n, i, j, jj, info, tmp, irej
      LOGICAL :: l_rejected

      find_unstable_mode = .FALSE.
      mode_eig = 0.0
      n = SIZE(dx)
      ALLOCATE(jacobian(n,n), eig_re(n), eig_im(n), eigvec(n,n), work(8*n), order(n))
      IF (.NOT.solve_linear_multi(dx_dx, dx_dg, jacobian)) RETURN
      CALL dgeev('N', 'V', n, jacobian, n, eig_re, eig_im, dummy, 1, eigvec, n, work, SIZE(work), info)
      IF (info /= 0) RETURN
      eig_abs_max = MAXVAL(SQRT(eig_re**2 + eig_im**2))
      IF (MINVAL(eig_re) >= -neg_eig_tolerance*eig_abs_max) RETURN

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
         ! dgeev stores Re(v) of a complex pair in the column of the member with Im(lambda) > 0
         IF (eig_im(j) < 0.0) CYCLE
         ! mode_dir = dX Re(v), mismatch = dG Re(v) - Re(lambda) dX Re(v)
         mode_dir = eigvec(1,j)*dx(1)
         mismatch = (-eigvec(1,j))*df(1)
         DO jj = 2, n
            mode_dir = mode_dir + eigvec(jj,j)*dx(jj)
            mismatch = mismatch - eigvec(jj,j)*df(jj)
         END DO
         dir_norm = SQRT(ABS(mode_dir.dot.mode_dir%apply_metric(l_dfpt)))
         IF (dir_norm == 0.0) CYCLE
         mismatch = mismatch - eig_re(j)*mode_dir
         mismatch_norm = SQRT(ABS(mismatch.dot.mismatch%apply_metric(l_dfpt)))
         IF (mismatch_norm > max_mode_inconsistency*ABS(eig_re(j))*dir_norm) CYCLE
         l_rejected = .FALSE.
         DO irej = 1, max_rejected_dirs
            IF (n_mix_steps - rejection_step(irej) > rejection_lifetime) CYCLE
            IF (ABS(w_rejected_dirs(irej).dot.mode_dir)/dir_norm > rejected_overlap) l_rejected = .TRUE.
         END DO
         IF (l_rejected) CYCLE
         mode_eig = eig_re(j)
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
