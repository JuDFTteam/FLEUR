!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_orbmom
   IMPLICIT NONE

   PRIVATE

   
   TYPE, PUBLIC :: t_orbmom
   !! Type representing orbital momentum coefficients.
   !! This type contains allocatable arrays for the orbital momentum components
   !! \(L_z\), \(L_+\), and \(L_-\), which are used in calculations involving
   !! radial functions and angular momentum quantum numbers.
      !> Coefficients for \(L_z\) (radfun, radfun, l, m).
      REAL, ALLOCATABLE :: lz(:, :, :, :)
      !> Coefficients for \(L_+\) (radfun, radfun, l, m).
      COMPLEX, ALLOCATABLE :: lp(:, :, :, :)
      !> Coefficients for \(L_-\) (radfun, radfun, l, m).
      COMPLEX, ALLOCATABLE :: lm(:, :, :, :)

   CONTAINS
      !> Initialize the orbital momentum coefficients.
      PROCEDURE, PASS :: init
      !> Calculate the orbital momentum contributions.
      PROCEDURE, PASS :: calc_orbmom
      !> Compute orbital hyperfine contribution from current k-point accumulation.
      PROCEDURE, PASS :: calc_hff_orbital
   END TYPE t_orbmom

CONTAINS

SUBROUTINE init(thisOrb, radfun, lmax)
   !! Initialize the orbital momentum coefficients.
   !! Allocates and initializes the arrays for \(L_z\), \(L_+\), and \(L_-\)
   !! based on the maximum radial function size and angular momentum quantum number.
      USE m_types_radfun
      IMPLICIT NONE

      !> The orbital momentum object to initialize.
      CLASS(t_orbmom), INTENT(INOUT) :: thisOrb
      !> Radial function information.
      TYPE(t_radfun), INTENT(IN) :: radfun
      !> Maximum angular momentum quantum number.
      INTEGER, INTENT(IN) :: lmax

      INTEGER :: dim1

      IF (ALLOCATED(thisOrb%lz)) DEALLOCATE(thisOrb%lz)
      IF (ALLOCATED(thisOrb%lp)) DEALLOCATE(thisOrb%lp)
      IF (ALLOCATED(thisOrb%lm)) DEALLOCATE(thisOrb%lm)

      dim1 = maxval(radfun%n_r)
      ALLOCATE(thisOrb%lz(dim1, dim1, 0:lmax, -lmax:lmax))
      ALLOCATE(thisOrb%lp(dim1, dim1, 0:lmax, -lmax:lmax))
      ALLOCATE(thisOrb%lm(dim1, dim1, 0:lmax, -lmax:lmax))
      thisOrb%lz = 0.0
      thisOrb%lp = 0.0
      thisOrb%lm = 0.0
   END SUBROUTINE init

   SUBROUTINE calc_orbmom(orb, abc, atoms, radfun, we, itype, jsp, clmom)
      !! Calculate the orbital momentum contributions.
      !! Computes the orbital momentum coefficients \(L_z\), \(L_+\), and \(L_-\)
      !! for a given atom type and spin channel. The results are accumulated into
      !! the provided `clmom` array.
      USE m_types_atoms
      USE m_types_abc
      USE m_types_radfun
      IMPLICIT NONE

      !> The orbital momentum object to calculate.
      CLASS(t_orbmom), INTENT(INOUT) :: orb
      !> Matching coefficients.
      TYPE(t_abc), INTENT(IN) :: abc
      !> Atomic structure information.
      TYPE(t_atoms), INTENT(IN) :: atoms
      !> Radial function information.
      TYPE(t_radfun), INTENT(IN) :: radfun
      !> Weights for occupied bands.
      REAL, INTENT(IN) :: we(:)
      !> Accumulated orbital momentum components (x, y, z).
      REAL, INTENT(INOUT) :: clmom(3)

      !> Atom type index.
      INTEGER, INTENT(IN) :: itype
      !> Spin channel index.
      INTEGER, INTENT(IN) :: jsp

      INTEGER :: natom, l, m, lm, i, j, jj
      REAL :: orbz, qmtlz(0:atoms%lmax(itype)), qmtlx(0:atoms%lmax(itype)), qmtly(0:atoms%lmax(itype))
      COMPLEX :: orbp, orbm

      CALL orb%init(radfun, atoms%lmax(itype))

      DO natom = 1, atoms%neq(itype)
         DO l = 0, atoms%lmax(itype)
            !     -----> sum over m
            DO m = -l, l
               lm = l*(l + 1) + m
                  !loop over radial functions
                  DO j = 1, radfun%n_r(l)
                     DO jj = 1, radfun%n_r(l)

                        ! coeff. for lz ->
                        orb%lz(j, jj, l, m) = orb%lz(j, jj, l, m) + sum(& ! sum over occupied bands
                                              we(:)*abc%cof(:, lm, j, natom)*CONJG(abc%cof(:, lm, jj, natom)))
                        ! coeff. for l+ <M'|l+|M>  
                        IF (m < l) THEN
                           orb%lp(j, jj, l, m) = orb%lp(j, jj, l, m) +sum( & ! sum over occupied bands
                                                 we(:)*abc%cof(:, lm, j, natom)*CONJG(abc%cof(:, lm + 1, jj, natom)))
                        ELSE
                           orb%lp(j, jj, l, m) = 0.0
                        END IF
                        ! coeff. for l- <M'|l-|M> 
                        IF (m > -l) THEN
                           orb%lm(j, jj, l, m) = orb%lm(j, jj, l, m) + sum(& ! sum over occupied bands
                                                 we(:)*abc%cof(:, lm, j, natom)*CONJG(abc%cof(:, lm - 1, jj, natom)))
                        ELSE
                           orb%lm(j, jj, l, m) = 0.0
                        END IF
                     end do
                  end do
               
            END DO
         END DO
      end do

      
      DO l = 0,atoms%lmax(itype)
         !--->    lm-decomposed density for each atom type
         orbz=0;orbp=0.0;orbm=0.0
         DO m = -l,l
            DO j=1,radfun%n_r(l)
               DO jj=1,radfun%n_r(l)
                  !lz
                  orbz = orbz+ m * orb%lz(j,jj,l,m) *radfun%integral(j,jj,l,jsp,jsp)
                  ! lx,ly
                  orbp = orbp+SQRT(REAL((l-m)*(l+m+1))) * ( orb%lp(j,jj,l,m) *radfun%integral(j,jj,l,jsp,jsp))
                  orbm = orbm+SQRT(REAL((l+m)*(l-m+1))) * ( orb%lm(j,jj,l,m) *radfun%integral(j,jj,l,jsp,jsp))
               enddo
            ENDDO
         enddo
         qmtlz(l) = orbz/ atoms%neq(itype)
         qmtlx(l) = 0.5*( REAL(orbp)+ REAL(orbm))/ atoms%neq(itype)
         qmtly(l) = 0.5*(AIMAG(orbp)-AIMAG(orbm))/ atoms%neq(itype)
      ENDDO
      clmom(1) = clmom(1)+ sum(qmtlx)
      clmom(2) = clmom(2)+ sum(qmtly)
      clmom(3) = clmom(3)+ sum(qmtlz)
   end subroutine

   SUBROUTINE calc_hff_orbital(orb, atoms, radfun, itype, jsp, hypFineOrbContribs)
      !! Compute the orbital hyperfine field contribution from the current k-point
      !! accumulation in orb%lz.  Must be called immediately after calc_orbmom
      !! (before the next call to calc_orbmom zeroes lz via orb%init).
      !! Both large and small radial components enter the 1/r³ matrix element.
      !! Results are accumulated (INOUT) so they sum over k-points.
      USE m_types_atoms
      USE m_types_radfun
      USE m_intgr, ONLY: intgr3
      IMPLICIT NONE

      CLASS(t_orbmom), INTENT(IN)    :: orb
      TYPE(t_atoms),   INTENT(IN)    :: atoms
      TYPE(t_radfun),  INTENT(IN)    :: radfun
      INTEGER,         INTENT(IN)    :: itype, jsp
      REAL,            INTENT(INOUT) :: hypFineOrbContribs(-1:3)

      INTEGER :: l, m, j, jj, lmax
      REAL    :: invR3Int, orbHff
      REAL    :: integrand(atoms%jmtd)

      lmax = min(3, atoms%lmax(itype))

      do l = 0, lmax
         orbHff = 0.0
         do m = -l, l
            do j = 1, radfun%n_r(l)
               do jj = 1, radfun%n_r(l)
                  integrand(1:atoms%jri(itype)) = &
                     (radfun%R(1:atoms%jri(itype), 1, j, l, jsp) * radfun%R(1:atoms%jri(itype), 1, jj, l, jsp) &
                    + radfun%R(1:atoms%jri(itype), 2, j, l, jsp) * radfun%R(1:atoms%jri(itype), 2, jj, l, jsp)) &
                    / atoms%rmsh(1:atoms%jri(itype), itype)**3
                  call intgr3(integrand(1:atoms%jri(itype)), atoms%rmsh(1:atoms%jri(itype), itype), &
                              atoms%dx(itype), atoms%jri(itype), invR3Int)
                  orbHff = orbHff + m * orb%lz(j, jj, l, m) * invR3Int
               end do
            end do
         end do
         orbHff = orbHff / atoms%neq(itype)
         hypFineOrbContribs(l) = hypFineOrbContribs(l) + orbHff
         hypFineOrbContribs(-1) = hypFineOrbContribs(-1) + orbHff
      end do
   END SUBROUTINE calc_hff_orbital

end module
