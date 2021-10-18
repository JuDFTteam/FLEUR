!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_usetup
   !-------------------------------------------------------------------+
   ! Sets up the quantities needed for the LDA+U subroutines:          |                                         |
   !     potential matrix: vs_mmp                                      |
   !     total energy contribution: e_ldau                             |
   !                                                  G.B. Oct. 2000   |
   !                                                                   |
   !     Extension to multiple U per atom type  G.M. 2017              |
   !-------------------------------------------------------------------+
   USE m_juDFT
   USE m_rotMMPmat
   USE m_dftUPotential
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE u_setup(atoms,input,noco,mpi,hub1data,inDen,pot,results)

      TYPE(t_mpi),      INTENT(IN)    :: mpi
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_hub1data), INTENT(IN)    :: hub1data
      TYPE(t_potden),   INTENT(IN)    :: inDen
      TYPE(t_potden),   INTENT(INOUT) :: pot
      TYPE(t_results),  INTENT(INOUT) :: results


      INTEGER :: itype,ispin,l,m,mp,i_u,n_u
      CHARACTER(len=2) :: l_type
      CHARACTER(len=9) :: l_form

      COMPLEX, ALLOCATABLE :: n_mmp(:,:,:,:)

      ! look, whether density matrix exists already:
      IF (ANY(ABS(inDen%mmpMat).GT.1e-12)) THEN

         n_u = atoms%n_u+atoms%n_hia

         !Rotate the density matrix if specified with phi or theta angles
         ALLOCATE(n_mmp,mold=inDen%mmpmat)
         DO i_u = 1, n_u
            n_mmp(:,:,i_u,:) = rotMMPmat(inDen%mmpmat(:,:,i_u,:),atoms%lda_u(i_u)%phi,&
                                         atoms%lda_u(i_u)%theta,0.0,atoms%lda_u(i_u)%l)
         ENDDO

         DO i_u = 1, n_u
            CALL dftUPotential(n_mmp(:,:,i_u,:), atoms, atoms%lda_u(i_u), input%jspins, &
                               any(noco%l_unrestrictMT).OR.input%ldauSpinoffd, &
                               pot%mmpMat(:,:,i_u,:))
         ENDDO

         IF (mpi%irank.EQ.0) THEN
            DO ispin = 1,SIZE(pot%mmpMat,4)
               WRITE (oUnit,'(a7,i3)') 'spin #',ispin
               DO i_u = 1, n_u
                  itype = atoms%lda_u(i_u)%atomType
                  l     = atoms%lda_u(i_u)%l

                  WRITE (l_type,'(i2)') 2*(2*l+1)
                  l_form = '('//l_type//'f12.7)'
                  WRITE (oUnit,'(a20,i3)') 'n-matrix for atom # ',itype
                  IF(i_u > atoms%n_u) WRITE(oUnit,'(A)') 'n-matrix calculated with DFT+Hubbard-1'
                  WRITE (oUnit,l_form) ((n_mmp(m,mp,i_u,ispin),m=-l,l),mp=-l,l)
                  WRITE (oUnit,'(a20,i3)') 'V-matrix for atom # ',itype
                  IF (atoms%lda_u(i_u)%l_amf) THEN
                     WRITE (oUnit,*) 'using the around-mean-field limit'
                  ELSE
                     WRITE (oUnit,*) 'using the atomic limit of LDA+U'
                  ENDIF
                  WRITE (oUnit,l_form) ((pot%mmpMat(m,mp,i_u,ispin),m=-l,l),mp=-l,l)
               END DO
            END DO
            WRITE (oUnit,*) results%e_ldau
         ENDIF
      ELSE
         IF (mpi%irank.EQ.0) THEN
            WRITE (*,*) 'no density matrix found ... skipping LDA+U'
         ENDIF
         pot%mmpMat = cmplx_0
         results%e_ldau = 0.0
      ENDIF

   END SUBROUTINE u_setup
END MODULE m_usetup
