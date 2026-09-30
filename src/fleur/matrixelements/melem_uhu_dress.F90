!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  A radial set dressed with the spherical Bessel function of one b, and the radial
!>  derivative of that product. It is what <u_{k+b1}|H|u_{k+b2}> needs on each side before
!>  the Hamiltonian of m_melem_uhu_radint can be applied, since a Bloch state at k+b carries
!>  j_L(|b| r) on top of the radial function.
!>
!>  Tabulated per modulus, for every radial slot, every l and every Bessel channel at once:
!>  the caller pairs two of these and never differentiates again.
!>
!>  THE BESSEL DERIVATIVE IS ANALYTIC, from the recurrence
!>
!>      j'_0(ar) = -a j_1(ar) ,   j'_L(ar) = a j_{L-1}(ar) - (L+1)/r j_L(ar)
!>
!>  on the j_L that sphbes has already returned. Differentiating the product numerically
!>  instead makes a four-point cubic stencil carry an oscillating factor, and on a
!>  logarithmic mesh that is where the accuracy goes. Only the radial function keeps a
!>  numerical derivative, and it is smooth.
!>
!>  THE FACTOR OF r, which is the trap. t_radfun stores r times the radial function, so what
!>  comes out here is r*u*j and r*d/dr(u*j) -- one power of r on the value and the same one
!>  power on the derivative, not two and not none. Written out, with s = r*u the stored one:
!>
!>      r d/dr(u j) = s' j - (s/r) j + s j'
!>
!>  The middle term is what a reader expects to be absent and is not: it comes from
!>  differentiating u = s/r rather than s. Dropping it costs a term of the same order as the
!>  one kept, and nothing downstream can see it.
!>
!>  SPLIT IN TWO ON PURPOSE. The arithmetic needs a mesh, one stored radial function and a
!>  modulus, and not one of those is a FLEUR type, so it is not handed one: t_radfun is
!>  unpacked in a thin layer above. What the split buys is that the half carrying the
!>  conventions depends on nothing, while t_radfun carries generate_radial_functions and
!>  through it reaches half the code base.
MODULE m_melem_uhu_dress
   USE m_juDFT
   USE m_types_atoms
   USE m_types_radfun
   USE m_sphbes
   USE m_differentiate, ONLY: difcub
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_uhu_dress
CONTAINS

   !> The thin layer: unpacks the types and hands each radial slot to the core.
   SUBROUTINE melem_uhu_dress(atoms, n, rk, rf, jspin, lwn, uj, duj)
      TYPE(t_atoms), INTENT(IN)  :: atoms
      INTEGER, INTENT(IN)        :: n         !> atom type
      REAL, INTENT(IN)           :: rk        !> |b|, the modulus this table belongs to
      TYPE(t_radfun), INTENT(IN) :: rf
      INTEGER, INTENT(IN)        :: jspin
      INTEGER, INTENT(IN)        :: lwn       !> lmax of this type
      REAL, INTENT(OUT) :: uj(:, :, :, 0:, 0:)   !> (jri,2,slot,l,L) = r*u*j
      REAL, INTENT(OUT) :: duj(:, :, :, 0:, 0:)  !> (jri,2,slot,l,L) = r*d/dr(u*j)

      INTEGER :: l, slot

      CALL timestart("melem_uhu_dress")
      uj = 0.0
      duj = 0.0
      DO l = 0, lwn
         DO slot = 1, rf%n_r(l)
            CALL melem_uhu_dress_core(atoms%jri(n), atoms%rmsh(1:, n), &
                                      rf%r(:, :, slot, l, jspin), rk, lwn, &
                                      uj(:, :, slot, l, 0:), duj(:, :, slot, l, 0:))
         END DO
      END DO
      CALL timestop("melem_uhu_dress")
   END SUBROUTINE melem_uhu_dress

   !> One radial function dressed with every Bessel channel of one modulus, and the radial
   !> derivative of each product. Plain arrays, no types: this is the part with the conventions
   !> in it, so it is the part that has to be checkable on its own.
   SUBROUTINE melem_uhu_dress_core(jri, rmsh, s, rk, lwn, uj, duj)
      INTEGER, INTENT(IN)  :: jri
      REAL,    INTENT(IN)  :: rmsh(:)      !> the radial mesh
      REAL,    INTENT(IN)  :: s(:, :)      !> (jri,2) the stored r*u, large and small
      REAL,    INTENT(IN)  :: rk           !> |b|
      INTEGER, INTENT(IN)  :: lwn
      REAL,    INTENT(OUT) :: uj(:, :, 0:)
      REAL,    INTENT(OUT) :: duj(:, :, 0:)

      REAL    :: jj(0:lwn, jri), djj(0:lwn, jri), jlpp(0:lwn), ds(jri, 2)
      REAL    :: xi
      INTEGER :: i, ll, c

      !> the Bessel table and its exact derivative, once for this modulus
      DO i = 1, jri
         xi = rmsh(i)
         CALL sphbes(lwn, rk*xi, jlpp)
         jj(0:lwn, i) = jlpp(0:lwn)
         djj(0, i) = -rk*jlpp(MIN(1, lwn))
         DO ll = 1, lwn
            djj(ll, i) = rk*jlpp(ll - 1) - REAL(ll + 1)/xi*jlpp(ll)
         END DO
      END DO

      CALL dstored(rmsh, s, jri, ds)
      DO ll = 0, lwn
         DO i = 1, jri
            xi = rmsh(i)
            DO c = 1, 2
               uj(i, c, ll) = s(i, c)*jj(ll, i)
               duj(i, c, ll) = ds(i, c)*jj(ll, i) - s(i, c)/xi*jj(ll, i) + s(i, c)*djj(ll, i)
            END DO
         END DO
      END DO
   END SUBROUTINE melem_uhu_dress_core

   !> d/dr of the stored r*u, both components, by a cubic through four mesh points. The three
   !> ends reuse the nearest interior stencil, which is what the mesh allows: a cubic needs
   !> four points and the first and last two have none to their side.
   SUBROUTINE dstored(rmsh, s, jri, ds)
      REAL, INTENT(IN)     :: rmsh(:)
      REAL, INTENT(IN)     :: s(:, :)     !> (jri,2) the stored r*u
      INTEGER, INTENT(IN)  :: jri
      REAL, INTENT(OUT)    :: ds(:, :)

      INTEGER :: i, c

      DO c = 1, 2
         ds(1, c) = difcub(rmsh(1:4), s(1:4, c), rmsh(1))
         DO i = 2, jri - 2
            ds(i, c) = difcub(rmsh(i - 1:i + 2), s(i - 1:i + 2, c), rmsh(i))
         END DO
         ds(jri - 1, c) = difcub(rmsh(jri - 3:jri), s(jri - 3:jri, c), rmsh(jri - 1))
         ds(jri, c) = difcub(rmsh(jri - 3:jri), s(jri - 3:jri, c), rmsh(jri))
      END DO
   END SUBROUTINE dstored

END MODULE m_melem_uhu_dress
