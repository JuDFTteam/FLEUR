MODULE m_checkMMPmat

   !Check whether the given density matrix makes sense (only diagonal)

   USE m_types
   USE m_juDFT
   USE m_constants
   use m_types_mixvector
   use m_mixing_history
   use m_mpi_bc_tool

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE checkMMPmat(indStart,indEnd,l_denmat_in_mixer,fmpi,atoms,input,outden,inden)

      INTEGER,             INTENT(IN)    :: indStart, indEnd
      LOGICAL,             INTENT(IN)    :: l_denmat_in_mixer
      TYPE(t_mpi),         INTENT(IN)    :: fmpi
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_potden),      INTENT(IN)    :: outden
      TYPE(t_potden),      INTENT(INOUT) :: inden

      REAL, PARAMETER :: lowBound = 0.0
      LOGICAL changed_elements
      INTEGER i_u,l,ispin,m
      REAL maxOcc

      changed_elements = .FALSE.

      maxOcc = (2 + 0.1)/input%jspins
      DO i_u = indStart, indEnd
         l = atoms%lda_u(i_u)%l
         !Check the diagonal elements
         DO ispin = 1, input%jspins
            DO m = -l,l
               IF(REAL(inden%mmpmat(m,m,i_u,ispin)).LT.lowBound) THEN
                  IF(REAL(outden%mmpmat(m,m,i_u,ispin)).GE.lowBound) THEN
                     inden%mmpmat(m,m,i_u,ispin) = 0.0
                     changed_elements = .TRUE.
                  ENDIF
               ELSE IF(REAL(inden%mmpmat(m,m,i_u,ispin)).GT.maxOcc) THEN
                  IF(REAL(outden%mmpmat(m,m,i_u,ispin)).LE.maxOcc) THEN
                     inden%mmpmat(m,m,i_u,ispin) = maxOcc
                     changed_elements = .TRUE.
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      CALL mpi_bc(changed_elements, 0, fmpi%mpi_comm)

      IF(changed_elements .AND. l_denmat_in_mixer) THEN
         IF(fmpi%irank.EQ.0) THEN
            WRITE(*,*) "Invalid elements in DFT+U matrix after mixing:"
         ENDIF
         CALL mixing_history_reset(fmpi)
         CALL mixvector_reset()
      ENDIF

      IF(.FALSE.) THEN
         WRITE(*,*) "-----------------------------------------------------------------"
         WRITE(*,*) "Using the Quasi-Newton methods for mixing and LDA+U"
         WRITE(*,*) "from the beginning of the SCF calculaion can be unstable."
         WRITE(*,*) "You can reset the mixing_history, use straight mixing for "
         WRITE(*,*) "the first iterations or use linear mixing for the density matrix"
         WRITE(*,*) "-----------------------------------------------------------------"
         CALL juDFT_error("Invalid elements in mmpmat", calledby="checkMMPmat")

      ENDIF

   END SUBROUTINE checkMMPmat

END MODULE m_checkMMPmat