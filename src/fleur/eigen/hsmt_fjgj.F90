!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_fjgj
   USE m_juDFT
   IMPLICIT NONE

   PRIVATE
   TYPE t_fjgj
      REAL, ALLOCATABLE    :: fj(:, :, :, :), gj(:, :, :, :)
   CONTAINS
      procedure :: alloc
      procedure :: calculate => hsmt_fjgj
      final :: finalize_fjgj
   END TYPE
   PUBLIC t_fjgj

CONTAINS
   subroutine alloc(fjgj, nvd, lmaxd, isp, noco)
      USE m_types
      CLASS(t_fjgj), INTENT(OUT) :: fjgj
      INTEGER, INTENT(IN)        :: nvd, lmaxd, isp
      TYPE(t_noco), INTENT(IN)   :: noco

      ALLOCATE (fjgj%fj(0:lmaxd, nvd, merge(1, isp, noco%l_noco):merge(2, isp, noco%l_noco), MERGE(2, 1, noco%l_ss)))
      ALLOCATE (fjgj%gj(0:lmaxd, nvd, merge(1, isp, noco%l_noco):merge(2, isp, noco%l_noco), MERGE(2, 1, noco%l_ss)))
      !$acc enter data copyin(fjgj)
      !$acc enter data create(fjgj%fj, fjgj%gj)
      !$acc kernels
      fjgj%fj = 0.0
      fjgj%gj = 0.0
      !$acc end kernels
   end subroutine

   subroutine finalize_fjgj(fjgj)
      TYPE(t_fjgj), INTENT(INOUT) :: fjgj
      !$acc exit data delete(fjgj%fj, fjgj%gj)
      !$acc exit data delete(fjgj)
   
      IF (ALLOCATED(fjgj%fj)) DEALLOCATE(fjgj%fj)
      IF (ALLOCATED(fjgj%gj)) DEALLOCATE(fjgj%gj)      
   end subroutine finalize_fjgj

   SUBROUTINE hsmt_fjgj(fjgj, input, atoms, cell, lapw, noco, usdus, n, ispin)
      !Calculate the fj&gj array which contain the part of the A,B matching coeff. depending on the
      !radial functions at the MT boundary as contained in usdus
      USE m_constants, ONLY: fpi_const
      USE m_sphbes
      USE m_types
      IMPLICIT NONE
      CLASS(t_fjgj), INTENT(INOUT) :: fjgj
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_cell), INTENT(IN)     :: cell
      TYPE(t_noco), INTENT(IN)     :: noco
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_lapw), INTENT(IN)     :: lapw
      TYPE(t_usdus), INTENT(IN)    :: usdus
      !     ..
      !     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: ispin, n

      !     ..
      !     .. Local Scalars ..

      INTEGER k, l, lo, intspin, jspin, jspinStart, jSpinEnd
      LOGICAL l_socfirst
      !     .. Local Arrays ..
      REAL ws(0:atoms%lmax(n),input%jspins)
      REAL, allocatable:: gs(:)
      REAL, allocatable:: gb(:, :), fb(:, :)
      LOGICAL apw(0:atoms%lmaxd)
      !     ..
      l_socfirst = noco%l_soc .AND. noco%l_noco .AND. (.NOT. noco%l_ss)

      DO l = 0, atoms%lmax(n)
         apw(l) = ANY(atoms%l_dulo(:atoms%nlo(n), n))
         IF ((input%l_useapw) .AND. (atoms%lapw_l(n) .GE. l)) apw(l) = .FALSE.
      END DO
      DO lo = 1, atoms%nlo(n)
         IF (atoms%l_dulo(lo, n)) apw(atoms%llo(lo, n)) = .TRUE.
      END DO

      IF (any(noco%l_constrained) .or. l_socfirst .OR. any(noco%l_unrestrictMT) .OR. any(noco%l_spinoffd_ldau)) THEN
         jspinStart = 1
         jspinEnd = input%jspins
      ELSE
         jspinEnd = ispin
         jspinStart = ispin
      END IF

      DO intspin = 1, MERGE(2, 1, noco%l_ss)
         gs = lapw%rk(:lapw%nv(intspin), intspin)*atoms%rmt(n)

         CALL d_sphbes(atoms%lmax(n), gs, fb, gb)
         !CPP_OMP PARALLEL DO default(none) shared(atoms, lapw, usdus, fjgj, cell,fb, gb, apw,  n, intspin, jspinStart, jspinEnd) &
         !CPP_OMP PRIVATE( l, jspin, ws)
         !$acc kernels copyin(lapw,lapw%nv,atoms,atoms%lmax,jspinStart,jspinEnd,cell,usdus)&
         !$acc copyin(usdus%us,usdus%usdus,usdus%us,usdus%duds,usdus%uds,usdus%dus,apw)&
         !$acc create(ws) present(fb,gb,apw,fjgj,fjgj%fj,fjgj%gj) default(none)
         !$acc loop private(k)
         DO k = 1, lapw%nv(intspin)
            !---> set up wronskians for the matching conditions for each ntype
            !$acc loop vector private(l, jspin)
            DO l = 0, atoms%lmax(n)
               DO jspin = jspinStart, jspinEnd
                  ws(l, jspin) = (fpi_const/SQRT(cell%omtil))/(usdus%uds(l, n, jspin)*usdus%dus(l, n, jspin) &
                                                               - usdus%us(l, n, jspin)*usdus%duds(l, n, jspin))
               END DO
            end do
            !$acc end loop
            DO jspin = jspinStart, jspinEnd
               !$acc loop private(l)
               DO l = 0, atoms%lmax(n)
                  IF (apw(l)) THEN
                     fjgj%fj(l, k, jspin, intspin) = 1.0*fpi_const/SQRT(cell%omtil)*fb(l, k)/usdus%us(l, n, jspin)
                     fjgj%gj(l, k, jspin, intspin) = 0.0
                  ELSE
                     fjgj%fj(l, k, jspin, intspin) = ws(l,jspin)*(usdus%uds(l, n, jspin)*lapw%rk(k, intspin)*gb(l, k) &
                                                                - usdus%duds(l, n, jspin)*fb(l, k))
                     fjgj%gj(l, k, jspin, intspin) = ws(l,jspin)*(usdus%dus(l, n, jspin)*fb(l, k) &
                                                                - usdus%us(l, n, jspin)*lapw%rk(k, intspin)*gb(l, k))
                  END IF
               END DO
               !$acc end loop
            END DO
         END DO
         !$acc end kernels
      END DO
      !$acc exit data delete(fb, gb)
      RETURN
   END SUBROUTINE hsmt_fjgj
END MODULE m_hsmt_fjgj
