!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This module generates the cmt coefficients and eigenvectors z   !
!     at all kpoints nkpt from the irreducible kpoints kpts%nkpt          !
!     and writes them out in cmt and z, respectively.                 !
!                                                 M.Betzinger(09/07)  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_gen_wavf

CONTAINS

   SUBROUTINE gen_wavf(kpts, sym, atoms, el_eig, ello_eig, cell, mpdata, vr0, &
                       hybdat, noco,nococonv, fmpi, input, jsp)

      ! nkpt       ::     number of all k-points
      USE m_types
      USE m_constants
      USE m_radfun
      USE m_radflo
      USE m_abcof
      USE m_trafo!, ONLY: waveftrafo_genwavf
      USE m_olap
      USE m_hyb_abcrot
      USE m_io_hybrid

      IMPLICIT NONE

      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      TYPE(t_mpdata), intent(in) :: mpdata
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_atoms), INTENT(IN)    :: atoms

      INTEGER, INTENT(IN)    :: jsp

      REAL, INTENT(IN)    :: vr0(:, :, :)!(jmtd,ntype,jspd)
      REAL, INTENT(IN)    :: el_eig(0:atoms%lmaxd, atoms%ntype)
      REAL, INTENT(IN)    :: ello_eig(:,:)

      ! local scalars
      INTEGER                 :: ilo
      INTEGER                 :: ikpt, itype, iop
      INTEGER                 :: l, ng

      INTEGER                 :: nodem, noded
      REAL                    :: wronk

      ! local arrays
      INTEGER                 :: rrot(3, 3, sym%nsym)
      INTEGER                 :: iarr(0:atoms%lmaxd, atoms%ntype)

      REAL                    :: vr(atoms%jmtd, atoms%ntype, input%jspins)
      REAL, ALLOCATABLE        :: u(:, :, :), du(:, :, :)

      REAL                    :: flo(atoms%jmtd, 2, atoms%nlod)
      REAL                    :: uuilon(atoms%nlod, atoms%ntype), duilon(atoms%nlod, atoms%ntype)
      REAL                    :: ulouilopn(atoms%nlod, atoms%nlod, atoms%ntype)

!     local arrays for abcof1
!      COMPLEX                 ::  a(nvd,0:lmd,natd,kpts%nkpt),b(nvd,0:lmd,natd,kpts%nkpt)

      TYPE(t_lapw)  :: lapw(kpts%nkptf)

      call timestart("gen_wavf")
      CALL hybdat%usdus%init(atoms, input%jspins)

      ! setup rotations in reciprocal space
      DO iop = 1, sym%nsym
         IF (iop <= sym%nop) THEN
            rrot(:, :, iop) = transpose(sym%mrot(:, :, sym%invtab(iop)))
         ELSE
            rrot(:, :, iop) = -rrot(:, :, iop - sym%nop)
         END IF
      END DO

      ! generate G-vectors, which fulfill |k+G|<rkmax
      ! for all k-points
      DO ikpt = 1, kpts%nkptf
         CALL lapw(ikpt)%init(input, noco,nococonv, kpts, atoms, sym, ikpt, cell)
      END DO

      ! set spherical component of the potential from the previous iteration vr
      vr = vr0

      ! calculate radial basis functions belonging to the
      ! potential vr stored in bas1 and bas2
      ! bas1 denotes the large component
      ! bas2    "     "  small component

      allocate(u(atoms%jmtd, 2, 0:atoms%lmaxd), &
                du(atoms%jmtd, 2, 0:atoms%lmaxd), &
                source=0.0)

      iarr = 2
      DO itype = 1, atoms%ntype
         IF (fmpi%irank == 0) WRITE (oUnit, FMT=8000) itype
         ng = atoms%jri(itype)
         DO l = 0, atoms%lmax(itype)
            CALL radfun(l, itype, jsp, el_eig(l, itype), vr(:, itype, jsp), &
                      atoms, u(:, :, l), du(:, :, l), hybdat%usdus, nodem, noded, wronk)
            IF (fmpi%irank == 0) WRITE (oUnit, FMT=8010) l, el_eig(l, itype), &
                               hybdat%usdus%us(l, itype, jsp), hybdat%usdus%dus(l, itype, jsp),&
                               nodem, hybdat%usdus%uds(l, itype, jsp), hybdat%usdus%duds(l, itype, jsp),&
                               noded, hybdat%usdus%ddn(l, itype, jsp), wronk

            hybdat%bas1(1:ng, 1, l, itype) = u(1:ng, 1, l)
            hybdat%bas2(1:ng, 1, l, itype) = u(1:ng, 2, l)
            hybdat%bas1(1:ng, 2, l, itype) = du(1:ng, 1, l)
            hybdat%bas2(1:ng, 2, l, itype) = du(1:ng, 2, l)

            hybdat%bas1_MT(1, l, itype) = hybdat%usdus%us(l, itype, jsp)
            hybdat%drbas1_MT(1, l, itype) = hybdat%usdus%dus(l, itype, jsp)
            hybdat%bas1_MT(2, l, itype) = hybdat%usdus%uds(l, itype, jsp)
            hybdat%drbas1_MT(2, l, itype) = hybdat%usdus%duds(l, itype, jsp)
         END DO

         IF (atoms%nlo(itype) >= 1) THEN
            CALL radflo(atoms, itype, jsp, ello_eig, vr(:, itype, jsp), u, du, fmpi, hybdat%usdus, uuilon, duilon, ulouilopn, flo)

            DO ilo = 1, atoms%nlo(itype)
               iarr(atoms%llo(ilo, itype), itype) = iarr(atoms%llo(ilo, itype), itype) + 1
               hybdat%bas1(1:ng, iarr(atoms%llo(ilo, itype), itype), atoms%llo(ilo, itype), itype) = flo(1:ng, 1, ilo)
               hybdat%bas2(1:ng, iarr(atoms%llo(ilo, itype), itype), atoms%llo(ilo, itype), itype) = flo(1:ng, 2, ilo)
               hybdat%bas1_MT(iarr(atoms%llo(ilo, itype), itype), atoms%llo(ilo, itype), itype) = hybdat%usdus%ulos(ilo, itype, jsp)
               hybdat%drbas1_MT(iarr(atoms%llo(ilo, itype), itype), atoms%llo(ilo, itype), itype) = hybdat%usdus%dulos(ilo, itype, jsp)
            END DO
         END IF
      END DO
      deallocate(u, du)

      ! consistency check
      IF (.not. all(iarr == mpdata%num_radfun_per_l)) call judft_error('gen_wavf: counting error')


8000  FORMAT(1x, /, /, ' wavefunction parameters for atom type', i3, ':', /, t32, 'radial function', t79, &
             'energy derivative', /, t3, 'l', t8, 'energy', t26, 'value', t39, 'derivative', t53, &
             'nodes', t68, 'value', t81, 'derivative', t95, 'nodes', t107, 'norm', t119, 'wronskian')
8010  FORMAT(i3, f10.5, 2(5x, 1p, 2e16.7, i5), 1p, 2e16.7)

      call timestop("gen_wavf")
   END SUBROUTINE gen_wavf
END MODULE m_gen_wavf
