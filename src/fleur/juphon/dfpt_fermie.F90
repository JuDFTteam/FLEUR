MODULE m_dfpt_fermie
   USE m_juDFT
#ifdef CPP_MPI
   USE mpi
#endif

CONTAINS
   SUBROUTINE dfpt_fermie(eig_id,dfpt_eig_id,fmpi,kpts,input,noco,results,results1)
      !! Calculate the perturbed occupation numbers from the unperturbed ones and the
      !! perturbed eigenenergies.
      !! This is only done for metals, i.e. systems where the smearing is not set
      !! to 0.
      !! Fermi-Dirac smearing is assumed.
      USE m_types
      USE m_constants
      USE m_eig66_io, ONLY : read_eig, write_eig

      IMPLICIT NONE

      TYPE(t_results), INTENT(INOUT) :: results, results1
      TYPE(t_mpi),     INTENT(IN)    :: fmpi
      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_noco),    INTENT(IN)    :: noco
      TYPE(t_kpts),    INTENT(IN)    :: kpts

      INTEGER, INTENT(IN) :: eig_id, dfpt_eig_id

      REAL    :: efermi, ef_num, ef_den, x
      INTEGER :: j, jsp, k, nspins, noccbd, noccbd_max

      REAL, ALLOCATABLE :: sxm(:,:,:)

#ifdef CPP_MPI
      INTEGER, PARAMETER :: comm = MPI_COMM_SELF
      INTEGER ierr
#endif

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

      RETURN
   END SUBROUTINE dfpt_fermie

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
