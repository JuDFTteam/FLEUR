MODULE m_hybrid_core

   ! read core radial wavefunctions from corebas
   ! corebas is written in cored.F
   ! (core basis functions can be read in once during an iteration)

CONTAINS
   SUBROUTINE corewf(atoms, jsp, input, dimension,&
  &                   vr, lmaxcd, maxindxc, mpi, lmaxc, nindxc, core1, core2, eig_c)
      USE m_types
      IMPLICIT NONE

      TYPE(t_mpi), INTENT(IN)   :: mpi
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_atoms), INTENT(IN)   :: atoms

      ! - scalars -
      INTEGER, INTENT(IN)        ::  jsp
      INTEGER, INTENT(IN)        :: lmaxcd
      INTEGER, INTENT(INOUT)     ::  maxindxc

      ! -arrays -
      INTEGER, INTENT(INOUT)     ::  lmaxc(:)
      INTEGER, INTENT(INOUT)     ::  nindxc(0:lmaxcd, atoms%ntype)

      REAL, INTENT(IN)       ::  vr(:, :, :)!(atoms%jmtd,atoms%ntypd,input%jspins)
      REAL, INTENT(INOUT)    ::  core1(:, :, 0:, :) !(atoms%jmtd,maxindxc,0:lmaxcd,atoms%ntype)
      REAL, INTENT(INOUT)     ::  core2(:, :, 0:, :) !(jmtd,maxindxc,0:lmaxcd,ntype)

      REAL, INTENT(INOUT)    ::  eig_c(maxindxc, 0:lmaxcd, atoms%ntype)

      ! - local scalars -
      INTEGER                   ::  ncstd
      INTEGER                   ::  ok, itype, i, j, idum, l

      REAL                      ::  rdum
      REAL                      ::  weight1, weight2
      ! - local arrays -
      INTEGER, ALLOCATABLE       ::  nindxcr(:, :)
      INTEGER, ALLOCATABLE       ::  l_qn(:, :), n_qn(:, :)

      REAL, ALLOCATABLE          ::  j_qn(:, :)
      REAL, ALLOCATABLE          ::  core1r(:, :, :, :), core2r(:, :, :, :)
      REAL, ALLOCATABLE          ::  core11r(:, :, :, :), core22r(:, :, :, :)
      REAL, ALLOCATABLE          ::  eig_cr(:, :, :)

      ncstd = maxval(atoms%ncst)
      ALLOCATE (nindxcr(0:ncstd, atoms%ntype), stat=ok)

      ! generate relativistic core wave functions( ->core1r,core2r )
      CALL calcorewf(dimension, input, jsp, atoms,&
     &                ncstd, vr,&
     &                lmaxc, nindxcr, core1r, core2r, eig_cr, mpi)

      nindxc = 0

      ! average over core states that only differ in j
      ! core functions with l-qn equal 0 doesnot change during the average process

      nindxc(0, :) = nindxcr(0, :)
      DO itype = 1, atoms%ntype
         core1(:, :nindxc(0, itype), 0, itype)&
      &    = core1r(:, 0, :nindxc(0, itype), itype)
         core2(:, :nindxc(0, itype), 0, itype)&
      &    = core2r(:, 0, :nindxc(0, itype), itype)
         eig_c(:nindxc(0, itype), 0, itype)&
      &    = eig_cr(0, :nindxc(0, itype), itype)
      END DO

      DO itype = 1, atoms%ntype
         DO l = 1, lmaxc(itype)
            weight1 = 2*(l - 0.5) + 1
            weight2 = 2*(l + 0.5) + 1
            IF (modulo(nindxcr(l, itype), 2) == 0) THEN
               DO i = 1, nindxcr(l, itype), 2
                  nindxc(l, itype) = nindxc(l, itype) + 1
                  core1(:atoms%jri(itype), nindxc(l, itype), l, itype) =&
         &                       (weight1*core1r(:atoms%jri(itype), l, i, itype) +&
         &                        weight2*core1r(:atoms%jri(itype), l, i + 1, itype))&
         &                       /(weight1 + weight2)
                  core2(:atoms%jri(itype), nindxc(l, itype), l, itype) =&
         &                       (weight1*core2r(:atoms%jri(itype), l, i, itype) + &
         &                        weight2*core2r(:atoms%jri(itype), l, i + 1, itype))&
         &                       /(weight1 + weight2)

                  eig_c(nindxc(l, itype), l, itype) =&
         &                       (weight1*eig_cr(l, i, itype) +&
         &                        weight2*eig_cr(l, i + 1, itype))&
         &                       /(weight1 + weight2)
               END DO
            ELSE
               DO i = 1, nindxcr(l, itype) - 1, 2
                  nindxc(l, itype) = nindxc(l, itype) + 1
                  core1(:atoms%jri(itype), nindxc(l, itype), l, itype) =&
         &                       (weight1*core1r(:atoms%jri(itype), l, i, itype) +&
         &                        weight2*core1r(:atoms%jri(itype), l, i + 1, itype))&
         &                       /(weight1 + weight2)
                  core2(:atoms%jri(itype), nindxc(l, itype), l, itype) =&
         &                       (weight1*core2r(:atoms%jri(itype), l, i, itype) + &
         &                        weight2*core2r(:atoms%jri(itype), l, i + 1, itype))&
         &                       /(weight1 + weight2)

                  eig_c(nindxc(l, itype), l, itype) =&
         &                       (weight1*eig_cr(l, i, itype) +&
         &                        weight2*eig_cr(l, i + 1, itype))&
         &                       /(weight1 + weight2)
               END DO
               nindxc(l, itype) = nindxc(l, itype) + 1
               core1(:atoms%jri(itype), nindxc(l, itype), l, itype)&
        &        = core1r(:atoms%jri(itype), l, nindxcr(l, itype), itype)
               core2(:atoms%jri(itype), nindxc(l, itype), l, itype)&
        &        = core2r(:atoms%jri(itype), l, nindxcr(l, itype), itype)
               eig_c(nindxc(l, itype), l, itype)&
        &        = eig_cr(l, nindxcr(l, itype), itype)
            END IF

         END DO
      END DO

      DEALLOCATE (nindxcr, core1r, core2r, eig_cr)

      IF (maxindxc /= maxval(nindxc))&
     &   STOP 'corewf: counting error nindxc'

   END SUBROUTINE corewf

   SUBROUTINE calcorewf(dimension, input, jspin, atoms,&
  &                     ncstd, vr,&
  &                     lmaxc, nindxcr, core1, core2, eig_c, mpi)

      USE m_intgr, ONLY: intgr3, intgr0, intgr1
      USE m_constants, ONLY: c_light
      USE m_setcor
      USE m_differ
      USE m_types
      IMPLICIT NONE

      TYPE(t_mpi), INTENT(IN)   :: mpi
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_atoms), INTENT(IN)   :: atoms

      !  - scalars -
      INTEGER, INTENT(IN) :: ncstd
      INTEGER, INTENT(IN) :: jspin
      INTEGER, INTENT(OUT):: lmaxc(:)

      !  - arrays -
      INTEGER, INTENT(OUT):: nindxcr(0:ncstd, atoms%ntype)
      REAL, INTENT(IN) :: vr(:, :, :)!(atoms%jmtd,atoms%ntypd,input%jspins)
      REAL, ALLOCATABLE :: core1(:, :, :, :), core2(:, :, :, :)
      REAL, ALLOCATABLE :: eig_c(:, :, :)

      !  - local scalars -
      INTEGER              :: i, j, itype, korb, ncmsh, nst, ierr
      REAL                 :: e, fj, fl, fn, t2, c, bmu, weight
      REAL                 :: d, dxx, rn, rnot, z, t1, rr
      LOGICAL, SAVE        :: first = .true.

      !  - local arrays -
      INTEGER              :: kappa(dimension%nstd), nprnc(dimension%nstd)
      REAL                 :: vrd(dimension%msh)
      REAL                 :: occ(dimension%nstd), occ_h(dimension%nstd, 2), a(dimension%msh), b(dimension%msh)
      REAL, ALLOCATABLE, SAVE:: vr0(:, :, :)

      !   - intrinsic functions -
      INTRINSIC exp, iabs, isign

      c = c_light(1.0)

      IF (first) THEN
         ALLOCATE (vr0(atoms%jmtd, atoms%ntype, input%jspins))
      END IF

      IF (input%frcor) THEN
         IF (first) THEN
            vr0 = vr
            first = .false.
         END IF
      ELSE
         vr0 = vr
         first = .false.
      END IF

      ! this loop determines the dimensions

      lmaxc = 0; nindxcr = 0
      DO itype = 1, atoms%ntype
         z = atoms%zatom(itype)
         dxx = atoms%dx(itype)
         bmu = 0.0
         CALL setcor(itype, input%jspins, atoms, input, bmu,&
      &               nst, kappa, nprnc, occ_h)

         IF ((bmu > 99.)) THEN
            occ(1:nst) = input%jspins*occ_h(1:nst, jspin)
         ELSE
            occ(1:nst) = occ_h(1:nst, 1)
         END IF
         rnot = atoms%rmsh(1, itype)
         d = exp(atoms%dx(itype))
         ncmsh = nint(log((atoms%rmt(itype) + 10.0)/rnot)/dxx + 1)
         ncmsh = min(ncmsh, dimension%msh)
         rn = rnot*(d**(ncmsh - 1))

         nst = atoms%ncst(itype)

         DO korb = 1, nst
            IF (occ(korb) > 0) THEN
               fn = nprnc(korb)
               fj = iabs(kappa(korb)) - .5e0
               weight = 2*fj + 1.e0
               IF (bmu > 99.) weight = occ(korb)
               fl = fj + (.5e0)*isign(1, kappa(korb))
               e = -2*(z/(fn + fl))**2

               nindxcr(NINT(fl), itype) = nindxcr(NINT(fl), itype) + 1
               lmaxc(itype) = max(lmaxc(itype), NINT(fl))
            ENDIF
         END DO

      END DO

      ALLOCATE (core1(atoms%jmtd, 0:maxval(lmaxc), maxval(nindxcr), atoms%ntype))
      ALLOCATE (core2(atoms%jmtd, 0:maxval(lmaxc), maxval(nindxcr), atoms%ntype))
      ALLOCATE (eig_c(0:maxval(lmaxc), maxval(nindxcr), atoms%ntype))

      core1 = 0; core2 = 0
      nindxcr = 0
      DO itype = 1, atoms%ntype
         z = atoms%zatom(itype)
         dxx = atoms%dx(itype)
         bmu = 0.0
         CALL setcor(itype, input%jspins, atoms, input, bmu, nst, kappa, nprnc, occ_h)

         IF ((bmu > 99.)) THEN
            occ(1:nst) = input%jspins*occ_h(1:nst, jspin)
         ELSE
            occ(1:nst) = occ_h(1:nst, 1)
         END IF
         rnot = atoms%rmsh(1, itype)
         d = exp(atoms%dx(itype))
         ncmsh = nint(log((atoms%rmt(itype) + 10.0)/rnot)/dxx + 1)
         ncmsh = min(ncmsh, dimension%msh)
         rn = rnot*(d**(ncmsh - 1))
         IF (mpi%irank == 0) THEN
            WRITE (6, FMT=8000) z, rnot, dxx, atoms%jri(itype)
         END IF
         DO j = 1, atoms%jri(itype)
            vrd(j) = vr0(j, itype, jspin)
         END DO

         if (input%l_core_confpot) THEN

            ! linear extension of the potential with slope t1 / a.u.
            t1 = 0.125
            t2 = vrd(atoms%jri(itype))/atoms%rmt(itype) - atoms%rmt(itype)*t1
            rr = atoms%rmt(itype)
         else
            t2 = vrd(atoms%jri(itype))/(atoms%jri(itype) - ncmsh)
         endif
         IF (atoms%jri(itype) < ncmsh) THEN
            DO i = atoms%jri(itype) + 1, ncmsh
               if (input%l_core_confpot) THEN
                  rr = d*rr
                  vrd(i) = rr*(t2 + rr*t1)
               else
                  vrd(i) = vrd(atoms%jri(itype)) + t2*(i - atoms%jri(itype))
               endif
!
            END DO
         END IF

         nst = atoms%ncst(itype)

         DO korb = 1, nst
            IF (occ(korb) > 0) THEN
               fn = nprnc(korb)
               fj = iabs(kappa(korb)) - .5e0
               weight = 2*fj + 1.e0
               IF (bmu > 99.) weight = occ(korb)
               fl = fj + (.5e0)*isign(1, kappa(korb))
               e = -2*(z/(fn + fl))**2
               CALL differ(fn, fl, fj, c, z, dxx, rnot, rn, d, ncmsh, vrd, e, a, b, ierr)

               nindxcr(NINT(fl), itype) = nindxcr(NINT(fl), itype) + 1

               core1(:atoms%jri(itype), NINT(fl), nindxcr(NINT(fl), itype), itype)&
        &                                                 = a(:atoms%jri(itype))
               core2(:atoms%jri(itype), NINT(fl), nindxcr(NINT(fl), itype), itype)&
        &                                                 = b(:atoms%jri(itype))

               eig_c(NINT(fl), nindxcr(NINT(fl), itype), itype) = e

               IF (mpi%irank == 0) THEN
                  WRITE (6, FMT=8010) fn, fl, fj, e, weight
               END IF
               IF (ierr /= 0) call judft_error('error in core-level routine')
            ENDIF
         END DO

      END DO

8000  FORMAT(/, /, 10x, 'z=', f4.0, 5x, 'r(1)=', e14.6, 5x, 'dx=', f8.6, 5x,&
      &       'm.t.index=', i4, /, 15x, 'n', 4x, 'l', 5x, 'j', 4x, 'energy', 7x,&
      &       'weight')
8010  FORMAT(12x, 2f5.0, f6.1, f10.4, f10.0)

   END SUBROUTINE calcorewf

   SUBROUTINE core_init(dimension, input, atoms, lmaxcd, maxindxc)

      USE m_intgr, ONLY: intgr3, intgr0, intgr1
      USE m_constants, ONLY: c_light
      USE m_setcor
      USE m_differ
      USE m_types
      IMPLICIT NONE

      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_input), INTENT(IN)       :: input
      TYPE(t_atoms), INTENT(IN)       :: atoms
      INTEGER, INTENT(OUT)            :: maxindxc, lmaxcd

      !  - local scalars -
      INTEGER              :: i, j, itype, korb, ncmsh, nst, ierr
      REAL                 :: e, fj, fl, fn, t, bmu, c
      REAL                 :: d, dxx, rn, rnot, z, t1, rr

      !  - local arrays -
      INTEGER              :: kappa(dimension%nstd), nprnc(dimension%nstd)
      INTEGER              :: nindxcr(0:dimension%nstd, atoms%ntype)
      REAL                 :: occ(dimension%nstd), occ_h(dimension%nstd, 2), a(dimension%msh), b(dimension%msh)
      INTEGER              :: lmaxc(atoms%ntype)

      !   - intrinsic functions -
      INTRINSIC exp, iabs, isign

      c = c_light(1.0)

      ! this loop determines the dimensions

      lmaxc = 0; nindxcr = 0
      DO itype = 1, atoms%ntype
         z = atoms%zatom(itype)
         dxx = atoms%dx(itype)
         bmu = 0.0
         CALL setcor(itype, input%jspins, atoms, input, bmu, nst, kappa, nprnc, occ_h)

         occ(1:nst) = occ_h(1:nst, 1)

         rnot = atoms%rmsh(1, itype)
         d = exp(atoms%dx(itype))
         ncmsh = nint(log((atoms%rmt(itype) + 10.0)/rnot)/dxx + 1)
         ncmsh = min(ncmsh, dimension%msh)
         rn = rnot*(d**(ncmsh - 1))

         nst = atoms%ncst(itype)

         DO korb = 1, nst
            IF (occ(korb) == 0) CYCLE
            fn = nprnc(korb)
            fj = iabs(kappa(korb)) - .5e0

            fl = fj + (.5e0)*isign(1, kappa(korb))
            e = -2*(z/(fn + fl))**2

            nindxcr(NINT(fl), itype) = nindxcr(NINT(fl), itype) + 1
            lmaxc(itype) = max(lmaxc(itype), NINT(fl))
         END DO

      END DO

      lmaxcd = maxval(lmaxc)
      maxindxc = maxval(nindxcr(0, :), nint((maxval(nindxcr)/2.0)))

   END SUBROUTINE core_init

END MODULE m_hybrid_core
