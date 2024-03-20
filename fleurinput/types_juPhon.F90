!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_juPhon
   USE m_judft
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_fleurinput_base) :: t_juPhon
      LOGICAL :: l_dfpt    = .FALSE. ! Phonon calculation on/off
      LOGICAL :: l_jpCheck = .FALSE. ! Check validity of input for a phonon run
      LOGICAL :: l_jpTest  = .FALSE. ! Run juPhon testset/inpu tests
      LOGICAL :: l_potout  = .FALSE. ! Write out potential
      LOGICAL :: l_eigout  = .FALSE. ! Write out eigenstuff
      LOGICAL :: l_symTsh  = .FALSE. ! Use symmetrized kinetic energy in SH
      LOGICAL :: l_symTdm  = .FALSE. ! /in dynamic matrix calculation
      LOGICAL :: l_bfkq    = .FALSE. ! Use backfolding of k+q to the first BZ
      INTEGER :: jplPlus   = 0       ! increase lmax for phonon vectors (lmax+1)
      REAL    :: kgqmax    = -1.0    ! Alternative maximum for |k+G+q|
      REAL    :: gqmax     = -1.0    ! Alternative maximum for |G+q|
      REAL    :: eps_pert  = 0.00001 ! Convergence criterion
      !REAL    :: eDiffcut  = 1e-12   ! Cutoff for energy differences
      !REAL    :: eDiffcut  = 1e-3   ! Cutoff for energy differences
      !REAL    :: eDiffcut  = 1e-7   ! Cutoff for energy differences
      REAL    :: eDiffcut  = 1e-5   ! Cutoff for energy differences
      REAL    :: fDiffcut  = 1e-7    ! Cutoff for occupation differences
      REAL    :: qpt_ph(3)           ! Debug q

      LOGICAL :: e1term = .TRUE.     ! Calculate the eigenenergy response
      LOGICAL :: f1term = .TRUE.     ! Calculate the occupation number response
      LOGICAL :: l_intp = .FALSE.    ! Interpolate the q-set onto another one
      LOGICAL :: l_band = .FALSE.    ! Interpolate the q-set to a bandstructure
      LOGICAL :: l_dos  = .FALSE.     ! Calculate the phonon density of states
      LOGICAL :: l_scf  = .TRUE.     ! Do a self-consistency run for dynmats
      INTEGER :: startq = 1          ! Start the q-loop at a specific point
      INTEGER :: stopq  = 0          ! Stop  the q-loop at a specific point
      INTEGER :: qmode  = 0          ! 0: Single-shot calculation for qlist
                                     ! 1: Reads q from fullsym_* input files

      REAL, ALLOCATABLE :: qvec(:,:)

      INTEGER :: singleQpt       = 1

      REAL    :: paPoX           = 1.0
      REAL    :: paPoY           = 1.0
      REAL    :: paPoZ           = 1.0
      LOGICAL :: harSw           = .TRUE.
      LOGICAL :: extSw           = .TRUE.
      LOGICAL :: xcSw            = .TRUE.
      LOGICAL :: calcEigenVec    = .TRUE.
      LOGICAL :: oneSternhCycle  = .FALSE.
      LOGICAL :: recipLengthUnit = .TRUE.
      LOGICAL :: onlyTests       = .FALSE.
      LOGICAL :: testsActivated  = .FALSE.
      INTEGER :: noPtsCon        = 30
      INTEGER :: numAddQs = 0

      LOGICAL :: testCompareGrVeff0FleurSw          = .FALSE.
      LOGICAL :: testVeff1Sw                        = .FALSE.
      LOGICAL :: testUnfoldStarsSw                  = .FALSE.
      LOGICAL :: testRadDerivativeSw                = .FALSE.
      LOGICAL :: testGauntCoeffSw                   = .FALSE.
      LOGICAL :: testGradLhExpandFuncSw             = .FALSE.
      LOGICAL :: testContGrVeff0Sw                  = .FALSE.
      LOGICAL :: testWarpingSw                      = .FALSE.
      LOGICAL :: testSternhHSMEtestSw               = .FALSE.
      LOGICAL :: testSternhSchroedPertTheoSw        = .FALSE.
      LOGICAL :: testz1Phi0ContSw                   = .FALSE.
      LOGICAL :: testRho1BasCorrSw                  = .FALSE.
      LOGICAL :: testPlotRho03Dsw                   = .FALSE.
      LOGICAL :: testRadSolSw                       = .FALSE.
      LOGICAL :: testKptsWeightSw                   = .FALSE.
      LOGICAL :: testCountValElecSw                 = .FALSE.
      LOGICAL :: testVeff0ContSw                    = .FALSE.
      LOGICAL :: testrho0ContSw                     = .FALSE.
      LOGICAL :: testBackRotMTCoordSysSw            = .FALSE.
      LOGICAL :: testRho1IRsw                       = .FALSE.
      LOGICAL :: testRho1MTsw                       = .FALSE.
      LOGICAL :: testPsi0ContSw                     = .FALSE.
      LOGICAL :: testOverlapSw                      = .FALSE.
      LOGICAL :: testGradRho0PathSw                 = .FALSE.
      LOGICAL :: testEii2LatPeriodQSw               = .FALSE.
      LOGICAL :: testVarphiHepsVarphiSw             = .FALSE.
      LOGICAL :: test1st2ndPulDynMatEps1            = .FALSE.
      LOGICAL :: test1st2ndPulDynMatCancel          = .FALSE.
      LOGICAL :: test3rdPulDynMatCancel             = .FALSE.
      LOGICAL :: testIntVeff1Rho1Val                = .FALSE.
      LOGICAL :: testGrPsiPsiMatElem                = .FALSE.
      LOGICAL :: testCompareSurfInt                 = .FALSE.
      LOGICAL :: testSplitMTSurfIntSterh            = .FALSE.
      LOGICAL :: testVeff1IRMESternh                = .FALSE.
      LOGICAL :: testEps1q0                         = .FALSE.
      LOGICAL :: testVeff1IRMatqBackFold            = .FALSE.
      LOGICAL :: testVeff1IRqLatPeriod              = .FALSE.
      LOGICAL :: testGrMatElemPsiHepsPsiGaussTheo   = .FALSE.
      LOGICAL :: testPsiHepsTildePsi                = .FALSE.
      LOGICAL :: testGoldsteinRemaining             = .FALSE.
      LOGICAL :: testR2orNotWfMtGradNgrNTensGrOvls  = .FALSE.
      LOGICAL :: testComp3ArgSFIntsSw               = .FALSE.
      LOGICAL :: testComp2ArgSFIntsSw               = .FALSE.
      LOGICAL :: testGoldsteinSurfSw                = .FALSE.
      LOGICAL :: testComp2ArgGrSFIntsSw             = .FALSE.
      LOGICAL :: testIRIntegralSw                   = .FALSE.
      LOGICAL :: testIR3rdMatElemSw                 = .FALSE.
      LOGICAL :: testActionHgrPhiSw                 = .FALSE.
      LOGICAL :: testXCintegrals                    = .FALSE.
      LOGICAL :: testEii2PsDens                     = .FALSE.

   CONTAINS
      PROCEDURE :: read_xml => read_xml_juPhon
      PROCEDURE :: mpi_bc => mpi_bc_juPhon
   END TYPE t_juPhon

   PUBLIC t_juPhon

CONTAINS

   SUBROUTINE mpi_bc_juPhon(this, mpi_comm, irank)
      USE m_mpi_bc_tool

      CLASS(t_juPhon), INTENT(INOUT) ::this
      INTEGER, INTENT(IN)           :: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL :: irank

      INTEGER :: rank

      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF

      CALL mpi_bc(this%l_potout, rank, mpi_comm)
      CALL mpi_bc(this%l_eigout, rank, mpi_comm)
      CALL mpi_bc(this%l_dfpt, rank, mpi_comm)
      CALL mpi_bc(this%e1term, rank, mpi_comm)
      CALL mpi_bc(this%f1term, rank, mpi_comm)
      CALL mpi_bc(this%l_intp, rank, mpi_comm)
      CALL mpi_bc(this%l_band, rank, mpi_comm)
      CALL mpi_bc(this%l_dos, rank, mpi_comm)
      CALL mpi_bc(this%l_scf, rank, mpi_comm)
      CALL mpi_bc(this%startq, rank, mpi_comm)
      CALL mpi_bc(this%stopq, rank, mpi_comm)
      CALL mpi_bc(this%qmode, rank, mpi_comm)
      CALL mpi_bc(this%singleQpt, rank, mpi_comm)
      CALL mpi_bc(this%qvec, rank, mpi_comm)

   END SUBROUTINE mpi_bc_juPhon

   SUBROUTINE read_xml_juPhon(this, xml)
      USE m_types_xml
      USE m_judft

      IMPLICIT NONE

      CLASS(t_juPhon), INTENT(INOUT) :: this
      TYPE(t_xml), INTENT(INOUT)     :: xml

      INTEGER::numberNodes
      CHARACTER(len=100) :: xPathA

      numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon')

      IF (numberNodes == 1) THEN
         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_dfpt')

         IF (numberNodes == 1) THEN
           this%l_dfpt    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_dfpt'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_jpCheck')

         IF (numberNodes == 1) THEN
           this%l_jpCheck = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_jpCheck'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_jpTest')

         IF (numberNodes == 1) THEN
           this%l_jpTest  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_jpTest'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_eigout')

         IF (numberNodes == 1) THEN
           this%l_potout  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_potout'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_potout')

         IF (numberNodes == 1) THEN
           this%l_eigout  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_eigout'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_symTsh')

         IF (numberNodes == 1) THEN
           this%l_symTsh  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_symTsh'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_symTdm')

         IF (numberNodes == 1) THEN
           this%l_symTdm  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_symTdm'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_bfkq')

         IF (numberNodes == 1) THEN
           this%l_bfkq    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_bfkq'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@jplPlus')

         IF (numberNodes == 1) THEN
           this%jplPlus    = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@jplPlus'))
         END IF

         IF ((this%jplPlus.NE.0).AND.(this%jplPlus.NE.1)) THEN
           CALL judft_error("Unreasonable increase of lmax!")
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@kgqmax')

         IF (numberNodes == 1) THEN
           this%kgqmax    = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@kgqmax'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@gqmax')

         IF (numberNodes == 1) THEN
           this%gqmax     = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@gqmax'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@eps_pert')

         IF (numberNodes == 1) THEN
           this%eps_pert  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@eps_pert'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@eDiffcut')

         IF (numberNodes == 1) THEN
           this%eDiffcut  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@eDiffcut'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@fDiffcut')

         IF (numberNodes == 1) THEN
           this%fDiffcut  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@fDiffcut'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@singleQpt')

         IF (numberNodes == 1) THEN
           this%singleQpt  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@singleQpt'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@e1term')

         IF (numberNodes == 1) THEN
           this%e1term  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@e1term'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@f1term')

         IF (numberNodes == 1) THEN
           this%f1term  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@f1term'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_intp')

         IF (numberNodes == 1) THEN
           this%l_intp  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_intp'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_band')

         IF (numberNodes == 1) THEN
           this%l_band  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_band'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_dos')

         IF (numberNodes == 1) THEN
           this%l_dos  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_dos'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_scf')

         IF (numberNodes == 1) THEN
           this%l_scf  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_scf'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@startq')

         IF (numberNodes == 1) THEN
           this%startq  = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@startq'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@stopq')

         IF (numberNodes == 1) THEN
           this%stopq  = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@stopq'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@qmode')

         IF (numberNodes == 1) THEN
           this%qmode  = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@qmode'))
         END IF

         this%qpt_ph(1) = 0.0
         this%qpt_ph(2) = 0.0
         this%qpt_ph(3) = 0.0

         allocate(this%qvec(0,0))
         this%qvec=xml%read_q_list('/fleurInput/output/juPhon/qVectors')
      ENDIF

   END SUBROUTINE read_xml_juPhon

END MODULE m_types_juPhon
