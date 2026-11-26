!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The radial integral of the spherical muffin-tin Hamiltonian between two products of a
!>  radial function and a spherical Bessel function,
!>
!>      < u j_1 | H_sph | u j_2 >
!>
!>  which is what turns a neighbour overlap into a <u_{k+b1}|H_k|u_{k+b2}>. The Bessel
!>  functions carry the two b vectors, so bra and ket live at different k points and neither
!>  is an eigenstate of H_k: the Hamiltonian has to be applied, not replaced by an eigenvalue.
!>
!>  Scalar relativistic, in the two-component form the radial functions already come in:
!>  Eq. (3.54) of the PhD thesis of P. Kurz. The mass enhancement is rebuilt per radial point
!>  from the energy parameter the radial functions were generated at.
!>
!>  H IS APPLIED AS AN OPERATOR AND THE OVERLAP IS TAKEN AFTERWARDS, which is what keeps the
!>  two sides honest. Writing the four large/small blocks out by hand invites a choice about
!>  how to spread H between bra and ket, and there is no choice to make: applying one
!>  operator to one side is unambiguous. Both orders are evaluated and averaged, so the result
!>  is Hermitian by construction rather than by a weight put in by hand.
!>
!>  Their difference is not noise. H is Hermitian over all space, but the sphere is cut at
!>  R_MT, so <a|Hb> - <Ha|b> is the surface term the truncation leaves behind. It comes out
!>  through 'asym' instead of being averaged away silently, because it measures how much the
!>  muffin-tin boundary costs and nothing else here can see it.
!>
!>  TWO CONVENTIONS THAT DO NOT ANNOUNCE THEMSELVES. vr is r times the potential, so it is
!>  divided by r here; and it is the spherical channel with its Y_00 factor already taken
!>  out, which the caller must not put back.
!>
!>  The limit it has to satisfy: at b1 = b2 = 0 the Bessel functions are one and their
!>  derivatives vanish, so uj reduces to u and the integral is the energy parameter times the
!>  norm t_radfun already tabulates, with asym down at the radial mesh error.
MODULE m_melem_uhu_radint
   USE m_constants, ONLY: c_light
   USE m_intgr, ONLY: intgr3
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_uhu_radint
CONTAINS

   !> One integral, for one pair of radial slots and one centrifugal channel. Called from the
   !> tabulation over lm pairs and neighbour pairs, so it allocates nothing: the integrands are
   !> automatic arrays and the caller owns everything else.
   SUBROUTINE melem_uhu_radint(jri, rmsh, dx, e, vr, uj1, uj2, duj1, duj2, lzb, &
                               integral, asym)
      INTEGER, INTENT(IN)  :: jri              !> radial points inside the sphere
      REAL,    INTENT(IN)  :: rmsh(:)          !> the radial mesh
      REAL,    INTENT(IN)  :: dx               !> its logarithmic increment
      REAL,    INTENT(IN)  :: e                !> energy parameter of the radial functions
      REAL,    INTENT(IN)  :: vr(:)            !> r * V(r), spherical channel, Y_00 taken out
      REAL,    INTENT(IN)  :: uj1(:, :)        !> (jri,2) u j at b1, large and small component
      REAL,    INTENT(IN)  :: uj2(:, :)        !> (jri,2) u j at b2
      REAL,    INTENT(IN)  :: duj1(:, :)       !> (jri,2) d/dr of uj1
      REAL,    INTENT(IN)  :: duj2(:, :)       !> (jri,2) d/dr of uj2
      INTEGER, INTENT(IN)  :: lzb              !> l of the centrifugal barrier
      REAL,    INTENT(OUT) :: integral         !> the Hermitian average
      REAL,    INTENT(OUT), OPTIONAL :: asym   !> <a|Hb> - <Ha|b>, the R_MT surface term

      REAL :: h1(jri, 2), h2(jri, 2), x(jri)
      REAL :: right, left

      CALL apply_hsra(jri, rmsh, e, vr, lzb, uj2, duj2, h2)
      x(:) = uj1(1:jri, 1)*h2(:, 1) + uj1(1:jri, 2)*h2(:, 2)
      CALL intgr3(x, rmsh, dx, jri, right)

      CALL apply_hsra(jri, rmsh, e, vr, lzb, uj1, duj1, h1)
      x(:) = h1(:, 1)*uj2(1:jri, 1) + h1(:, 2)*uj2(1:jri, 2)
      CALL intgr3(x, rmsh, dx, jri, left)

      integral = 0.5*(right + left)
      IF (PRESENT(asym)) asym = right - left
   END SUBROUTINE melem_uhu_radint

   !> The spherical scalar-relativistic Hamiltonian acting on one two-component radial
   !> function. The whole operator lives here and nowhere else, which is why the caller has
   !> no weights to choose.
   SUBROUTINE apply_hsra(jri, rmsh, e, vr, lzb, uj, duj, huj)
      INTEGER, INTENT(IN)  :: jri
      REAL,    INTENT(IN)  :: rmsh(:), vr(:), e
      INTEGER, INTENT(IN)  :: lzb
      REAL,    INTENT(IN)  :: uj(:, :), duj(:, :)
      REAL,    INTENT(OUT) :: huj(:, :)

      REAL    :: c, c2, cin2, ll, xi, vv, mm
      INTEGER :: i

      c = c_light(1.0)
      c2 = c*c
      cin2 = 1.0/c2
      ll = REAL(lzb*(lzb + 1))

      DO i = 1, jri
         xi = rmsh(i)
         vv = vr(i)/xi
         !> the scalar-relativistic mass enhancement, at the energy the functions were made at
         mm = 1.0 + 0.5*cin2*(e - vv)
         huj(i, 1) = uj(i, 1)*(0.5/mm*ll/xi/xi + vv) - c*(2.0*uj(i, 2)/xi + duj(i, 2))
         huj(i, 2) = uj(i, 2)*(-2.0*c2 + vv) + c*duj(i, 1)
      END DO
   END SUBROUTINE apply_hsra

END MODULE m_melem_uhu_radint
