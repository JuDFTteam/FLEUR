MODULE m_dfpt_fermie
   USE m_juDFT
#ifdef CPP_MPI
   USE mpi
#endif

CONTAINS
   SUBROUTINE dfpt_fermie(eig_id,dfpt_eig_id,fmpi,kpts,input,noco,results,results1)

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
      INTEGER :: j, jsp, k, nspins, noccbd

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

      DO jsp = 1, nspins
         DO k = 1, kpts%nkpt
            IF (fmpi%irank == 0) THEN
               IF (input%eig66(1)) CALL read_eig(dfpt_eig_id,k,jsp,neig=results1%neig(k,jsp),eig=results1%eig(:,k,jsp))
            END IF
#ifdef CPP_MPI
            CALL MPI_BARRIER(fmpi%mpi_comm,ierr)
#endif
         END DO
      END DO

      IF (fmpi%irank == 0) THEN
         efermi = results%ef
         results1%ef = 0.0
         ef_num = 0.0
         ef_den = 0.0

         DO jsp = 1, nspins
            DO k = 1, kpts%nkpt
               noccbd  = COUNT(results%w_iks(:,k,jsp)*2.0/input%jspins>1.e-8)
               DO j = 1, noccbd
                  !x = (results%eig(j,k,jsp)-efermi)/input%tkb
                  !sxm(j,k,jsp) = sfermi(-x)
                  !ef_num = ef_num + results%w_iks(j,k,jsp) * sxm(j,k,jsp) * results1%eig(j,k,jsp)
                  !ef_den = ef_den + results%w_iks(j,k,jsp) * sxm(j,k,jsp)
                  x = results%w_iks(j,k,jsp) * (1 - results%w_iks(j,k,jsp)/kpts%wtkpt(k))
                  sxm(j,k,jsp) = x
                  ef_num = ef_num + x * results1%eig(j,k,jsp)
                  ef_den = ef_den + x
               END DO
            END DO
         END DO

         IF (ABS(ef_den)>1e-12) THEN
            results1%ef = ef_num/ef_den
         ELSE
            results1%ef = 0.0
         END IF

         !results1%w_iks(:noccbd,:,1:nspins) = -results%w_iks(:noccbd,:,1:nspins) &
         !                             * sxm(:noccbd,:,1:nspins) &
         !                             * (results1%eig(:noccbd,:,1:nspins)-results1%ef)/input%tkb
         results1%w_iks(:noccbd,:,1:nspins) = - sxm(:noccbd,:,1:nspins) &
                                      * (results1%eig(:noccbd,:,1:nspins)-results1%ef)/input%tkb
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
