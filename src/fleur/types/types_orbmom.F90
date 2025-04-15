!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_orbmom
   IMPLICIT NONE

   PRIVATE

   TYPE,PUBLIC:: t_orbmom
      REAL, ALLOCATABLE    :: lz(:, :, :, :) !!(radfun,radfun,l,m)
      COMPLEX, ALLOCATABLE :: lp(:, :, :, :) !!(radfun,radfun,l,m)
      COMPLEX, ALLOCATABLE :: lm(:, :, :, :) !!(radfun,radfun,l,m)

   CONTAINS
      PROCEDURE, PASS :: calc_orbmom,init
   END TYPE t_orbmom

CONTAINS

   SUBROUTINE init(thisOrb, radfun, lmax)

      USE m_types_radfun
      IMPLICIT NONE

      CLASS(t_orbmom), INTENT(INOUT)    :: thisOrb
      TYPE(t_radfun), INTENT(IN)     :: radfun
      INTEGER, INTENT(IN)            :: lmax

      INTEGER                        :: dim1

      IF (ALLOCATED(thisOrb%lz)) DEALLOCATE (thisOrb%lz)
      IF (ALLOCATED(thisOrb%lp)) DEALLOCATE (thisOrb%lp)
      IF (ALLOCATED(thisOrb%lm)) DEALLOCATE (thisOrb%lm)

      dim1 = maxval(radfun%n_r)
      allocate (thisOrb%lz(dim1, dim1, 0:lmax, -lmax:lmax))
      allocate (thisOrb%lp(dim1, dim1, 0:lmax, -lmax:lmax))
      allocate (thisOrb%lm(dim1, dim1, 0:lmax, -lmax:lmax))
      thisOrb%lz = 0.0
      thisOrb%lp = 0.0
      thisOrb%lm = 0.0

   END SUBROUTINE init

   subroutine calc_orbmom(orb, abc, atoms, radfun,we, itype,jsp,clmom)
      use m_types_atoms
      use m_types_abc
      use m_types_radfun
      class(t_orbmom), intent(inout):: orb
      type(t_abc), intent(in):: abc
      type(t_atoms), intent(in):: atoms
      type(t_radfun),intent(in):: radfun
      real, intent(in):: we(:)
      real,intent(out)::clmom(3)

      integer, intent(in) :: itype,jsp
      integer:: natom, l, m, lm, i, j, jj
      
      real:: orbz,qmtlz(0:atoms%lmax(itype)),qmtlx(0:atoms%lmax(itype)),qmtly(0:atoms%lmax(itype))
      complex:: orbp,orbm
      
      call orb%init(radfun,atoms%lmax(itype))
      
      DO natom = 1, atoms%neq(itype)
         DO l = 0, atoms%lmax(itype)
            !     -----> sum over m
            DO m = -l, l
               lm = l*(l + 1) + m
               !     -----> sum over occupied bands
               DO i = 1, size(we)
                  DO j = 1, radfun%n_r(l)
                     DO jj = 1, radfun%n_r(l)

                        ! coeff. for lz ->
                        orb%lz(j, jj, l, m) = orb%lz(j, jj, l, m) + &
                                              we(i)*abc%cof(i, lm, j, natom)*CONJG(abc%cof(i, lm, jj, natom))
                        ! coeff. for l+ <M'|l+|M> with respect to M ->
                        IF (m .NE. l) THEN
                           orb%lp(j, jj, l, m) = orb%lp(j, jj, l, m) + &
                                                 we(i)*abc%cof(i, lm, j, natom)*CONJG(abc%cof(i, lm + 1, jj, natom))
                        ELSE
                           orb%lp(j, jj, l, m) = 0.0
                        END IF
                        ! coeff. for l- <M'|l-|M> with respect to M ->
                        IF (m .NE. -l) THEN
                           orb%lm(j, jj, l, m) = orb%lm(j, jj, l, m) + &
                                                 we(i)*abc%cof(i, lm, j, natom)*CONJG(abc%cof(i, lm - 1, jj, natom))
                        ELSE
                           orb%lm(j, jj, l, m) = 0.0
                        END IF
                     end do
                  end do
               END DO
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
                  orbm = orbp+SQRT(REAL((l+m)*(l-m+1))) * ( orb%lm(j,jj,l,m) *radfun%integral(j,jj,l,jsp,jsp))
               enddo
            ENDDO
         ENDDO
         qmtlz(l)  = orbz/ atoms%neq(itype)
         qmtlx(l) = 0.5*( REAL(orbp)+ REAL(orbm))/ atoms%neq(itype)
         qmtly(l) = 0.5*(AIMAG(orbp)-AIMAG(orbm))/ atoms%neq(itype)
      ENDDO

      clmom(1) = sum(qmtlx)
      clmom(2) = sum(qmtly)
      clmom(3) = sum(qmtlz)
   end subroutine
end module
