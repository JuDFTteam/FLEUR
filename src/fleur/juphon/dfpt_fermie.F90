!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_fermie
   USE m_juDFT
#ifdef CPP_MPI
   USE mpi
#endif

CONTAINS
   SUBROUTINE dfpt_fermie(fmpi,kpts,input,noco,results,results1)
      !! Calculate the perturbed occupation numbers from the unperturbed ones and the
      !! perturbed eigenenergies.
      !! This is only done for metals, i.e. systems with partially occupied bands.
      !!
      !! The occupation-number derivative needed here is fundamentally different
      !! depending on how the ground state's Brillouin-zone integration was done:
      !! Fermi-Dirac smearing (hist/gauss/tria, all sharing input%tkb) gives a
      !! diagonal derivative -- the occupation of (n,k) depends only on its own
      !! eigenvalue. Tetrahedron integration (tetra) has no smearing width and
      !! is not diagonal -- the weight of one tetrahedron corner depends on all
      !! of that tetrahedron's corner eigenvalues. Dispatch to the matching
      !! implementation accordingly.
      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_results), INTENT(INOUT) :: results, results1
      TYPE(t_mpi),     INTENT(IN)    :: fmpi
      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_noco),    INTENT(IN)    :: noco
      TYPE(t_kpts),    INTENT(IN)    :: kpts

      IF (input%bz_integration==BZINT_METHOD_TETRA) THEN
         CALL dfpt_fermie_tetra(fmpi,kpts,input,noco,results,results1)
      ELSE
         CALL dfpt_fermie_hist(fmpi,kpts,input,noco,results,results1)
      END IF

   END SUBROUTINE dfpt_fermie

   SUBROUTINE dfpt_fermie_hist(fmpi,kpts,input,noco,results,results1)
      !! Fermi-Dirac smearing derivative of the perturbed occupation numbers,
      !! for bz_integration==hist/gauss/tria (all use input%tkb-based smearing
      !! for the ground-state occupations, so the same derivative applies).
      USE m_types

      IMPLICIT NONE

      TYPE(t_results), INTENT(INOUT) :: results, results1
      TYPE(t_mpi),     INTENT(IN)    :: fmpi
      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_noco),    INTENT(IN)    :: noco
      TYPE(t_kpts),    INTENT(IN)    :: kpts

      REAL    :: efermi, ef_num, ef_den, x
      INTEGER :: j, jsp, k, nspins, noccbd, noccbd_max

      REAL, ALLOCATABLE :: sxm(:,:,:)

      IF (noco%l_noco) THEN
         nspins = 1
      ELSE
         nspins = input%jspins
      END IF

      ALLOCATE(sxm(MAXVAL(results%neig),kpts%nkpt,nspins))

      IF (fmpi%irank == 0) THEN
         efermi = results%ef
         results1%ef = 0.0
         ef_num = 0.0
         ef_den = 0.0
         noccbd_max = 0
         sxm = 0.0

         DO jsp = 1, nspins
            DO k = 1, kpts%nkpt
               noccbd  = COUNT(results%w_iks(:,k,jsp)*2.0/input%jspins>1.e-8)
               IF (noccbd > noccbd_max ) noccbd_max = noccbd
               DO j = 1, noccbd
                  x = (results%eig(j,k,jsp)-efermi)/input%tkb
                  sxm(j,k,jsp) = sfermi(-x)
                  ef_num = ef_num + results%w_iks(j,k,jsp) * sxm(j,k,jsp) * results1%eig(j,k,jsp)
                  ef_den = ef_den + results%w_iks(j,k,jsp) * sxm(j,k,jsp)
               END DO
            END DO
         END DO

         IF (ABS(ef_den)>1e-12) THEN
            results1%ef = ef_num/ef_den
         ELSE
            results1%ef = 0.0
         END IF

         results1%w_iks(:noccbd_max,:,1:nspins) = -results%w_iks(:noccbd_max,:,1:nspins) &
                                            * sxm(:noccbd_max,:,1:nspins) &
                                            * (results1%eig(:noccbd_max,:,1:nspins)-results1%ef)/input%tkb
      END IF

   END SUBROUTINE dfpt_fermie_hist

   SUBROUTINE dfpt_fermie_tetra(fmpi,kpts,input,noco,results,results1)
      !! Tetrahedron-method equivalent of dfpt_fermie_hist.
      !!
      !! The tetrahedron weight assigned to one corner of a tetrahedron depends
      !! on all corner eigenvalues of that tetrahedron (not just its own), so
      !! the full corner-to-corner Jacobian of the tetrahedron weight is needed,
      !! not just a diagonal (own-eigenvalue) derivative. This routine builds
      !! that Jacobian by finite-differencing the analytic per-tetrahedron
      !! weight (tetraWeight, plus the Bloechl correction where enabled) corner
      !! by corner, and combines it with the already-known perturbed
      !! eigenvalues (results1%eig) at all corners of every tetrahedron to
      !! solve for the perturbed Fermi energy (from charge neutrality) and the
      !! perturbed occupation-number weights.
      USE m_types

      IMPLICIT NONE

      TYPE(t_results), INTENT(INOUT) :: results, results1
      TYPE(t_mpi),     INTENT(IN)    :: fmpi
      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_noco),    INTENT(IN)    :: noco
      TYPE(t_kpts),    INTENT(IN)    :: kpts

      REAL, PARAMETER :: del = 1.0e-4

      REAL    :: ef_num, ef_den, efermi
      INTEGER :: jsp, nspins, neig, noccbd_max, nc
      INTEGER :: itet, icorn, j, iband, ikptTarget
      REAL    :: vol(kpts%ntet)
      REAL    :: etetra(SIZE(kpts%ntetra,1))
      REAL    :: etetraP(SIZE(kpts%ntetra,1)), etetraM(SIZE(kpts%ntetra,1))
      REAL    :: dwdeps, sumdwdeps, crossterm

      REAL, ALLOCATABLE :: crossterm_acc(:,:,:), sumdwdeps_acc(:,:,:)

      IF (noco%l_noco) THEN
         nspins = 1
      ELSE
         nspins = input%jspins
      END IF

      IF (fmpi%irank == 0) THEN
         nc  = SIZE(kpts%ntetra,1)
         vol = kpts%voltet(:)/kpts%ntet

         ALLOCATE(crossterm_acc(MAXVAL(results%neig),kpts%nkpt,nspins))
         ALLOCATE(sumdwdeps_acc(MAXVAL(results%neig),kpts%nkpt,nspins))
         crossterm_acc = 0.0
         sumdwdeps_acc = 0.0

         DO jsp = 1, nspins
            neig   = MINVAL(results%neig(:,jsp))
            efermi = results%ef

            !$OMP parallel do default(none) &
            !$OMP shared(kpts,nc,neig,vol,input,results,results1,jsp,efermi) &
            !$OMP shared(crossterm_acc,sumdwdeps_acc) &
            !$OMP private(itet,iband,etetra,icorn,ikptTarget,sumdwdeps,crossterm) &
            !$OMP private(j,etetraP,etetraM,dwdeps) &
            !$OMP collapse(2)
            DO itet = 1, kpts%ntet
               DO iband = 1, neig

                  etetra = results%eig(iband,kpts%ntetra(:,itet),jsp)
                  IF( ALL(etetra>efermi+1e-8) .AND. .NOT.input%l_bloechl ) CYCLE

                  DO icorn = 1, nc
                     ikptTarget = kpts%ntetra(icorn,itet)

                     sumdwdeps = 0.0
                     crossterm = 0.0
                     DO j = 1, nc
                        etetraP = etetra
                        etetraP(j) = etetra(j) + del
                        etetraM = etetra
                        etetraM(j) = etetra(j) - del

                        dwdeps = ( cornerWeight(efermi,etetraP,icorn,vol(itet),input%film,input%l_bloechl) &
                                 - cornerWeight(efermi,etetraM,icorn,vol(itet),input%film,input%l_bloechl) ) / (2.0*del)

                        sumdwdeps = sumdwdeps + dwdeps
                        crossterm = crossterm + dwdeps * results1%eig(iband,kpts%ntetra(j,itet),jsp)
                     END DO

                     !$OMP critical
                     crossterm_acc(iband,ikptTarget,jsp) = crossterm_acc(iband,ikptTarget,jsp) + crossterm
                     sumdwdeps_acc(iband,ikptTarget,jsp) = sumdwdeps_acc(iband,ikptTarget,jsp) + sumdwdeps
                     !$OMP end critical
                  END DO
               END DO
            END DO
            !$OMP end parallel do
         END DO

         ef_num = 0.0
         ef_den = 0.0
         noccbd_max = 0
         DO jsp = 1, nspins
            neig = MINVAL(results%neig(:,jsp))
            IF (neig > noccbd_max) noccbd_max = neig
            ef_num = ef_num + SUM(crossterm_acc(:neig,:,jsp))
            ef_den = ef_den + SUM(sumdwdeps_acc(:neig,:,jsp))
         END DO

         IF (ABS(ef_den)>1e-12) THEN
            results1%ef = ef_num/ef_den
         ELSE
            results1%ef = 0.0
         END IF

         results1%w_iks = 0.0
         results1%w_iks(:noccbd_max,:,1:nspins) = crossterm_acc(:noccbd_max,:,1:nspins) &
                                       - sumdwdeps_acc(:noccbd_max,:,1:nspins) * results1%ef
      END IF

   CONTAINS

      PURE REAL FUNCTION cornerWeight(efermi,etetra,icorn,vol,film,l_bloechl) RESULT(w)
         !! Tetrahedron weight (+ Bloechl correction, if enabled) assigned to
         !! corner icorn, for the (unsorted, corner-index-ordered) eigenvalues
         !! etetra at a fixed Fermi energy efermi.
         USE m_tetsrt
         USE m_tetraWeight
         USE m_bloechl

         REAL,    INTENT(IN) :: efermi
         REAL,    INTENT(IN) :: etetra(:)
         INTEGER, INTENT(IN) :: icorn
         REAL,    INTENT(IN) :: vol
         LOGICAL, INTENT(IN) :: film, l_bloechl

         INTEGER :: ind(SIZE(etetra)), icornSorted, i

         ind = tetsrt(etetra)
         icornSorted = 1
         DO i = 1, SIZE(etetra)
            IF(ind(i)==icorn) icornSorted = i
         END DO

         w = tetraWeight(efermi,etetra(ind),icornSorted,vol,film)
         IF(l_bloechl) w = w + bloechl(efermi,etetra(ind),icornSorted,vol,film)

      END FUNCTION cornerWeight

   END SUBROUTINE dfpt_fermie_tetra

   REAL FUNCTION sfermi(x)
      !! Returns the Fermi-Dirac function
      !! $$s(x)=(e^{x}+1)^{-1}$$
      !! for \(x=(\epsilon_{\nu\boldsymbol{k}}-E_{F})/(k_{B}T)\).

      REAL, INTENT(IN) :: x

      REAL :: expo

      expo = EXP(x)

      sfermi = 1.0/(expo+1.0)

      RETURN

   END FUNCTION sfermi

END MODULE m_dfpt_fermie
