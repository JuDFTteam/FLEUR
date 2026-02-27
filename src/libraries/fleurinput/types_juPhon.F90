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
      !LOGICAL :: l_jpTest  = .FALSE. ! Run juPhon testset/inpu tests
      !LOGICAL :: l_potout  = .FALSE. ! Write out potential
      !LOGICAL :: l_eigout  = .FALSE. ! Write out eigenstuff
      REAL    :: eDiffcut  = 1e-5   ! Cutoff for energy differences
      REAL    :: fDiffcut  = 1e-7    ! Cutoff for occupation differences
      REAL    :: qlim      = 1./100     ! qlim value
      REAL    :: gmaxzLocal   = 0.0  ! Local Gmaxz cutoff for film 

      LOGICAL :: l_intp = .FALSE.    ! Interpolate the q-set onto another one
      LOGICAL :: l_band = .FALSE.    ! Interpolate the q-set to a bandstructure
      LOGICAL :: l_dos  = .FALSE.     ! Calculate the phonon density of states
      LOGICAL :: l_scf  = .TRUE.     ! Do a self-consistency run for dynmats
      LOGICAL :: l_elph = .FALSE.    ! Calculate electron-phonon matrix elements
      LOGICAL :: l_sumrule  = .FALSE. ! Apply sumrule for dynmats
      LOGICAL :: l_rm_qhdf  = .TRUE. ! Remove q*hdf files, after convergence
      INTEGER :: startq = 1          ! Start the q-loop at a specific point
      INTEGER :: stopq  = 0          ! Stop  the q-loop at a specific point
      INTEGER :: qmode  = 0          ! 0: Single-shot calculation for qlist
                                     ! 1: Reads q from fullsym_* input files
      INTEGER :: i_integration = 1   ! choose integration scheme for ph-linewidth (experimental)
      REAL    :: smearingGauss = 1e-7 ! Gaussian smearing for pot-response in Film DFPT 
      LOGICAL :: l_phonon = .TRUE.    
      LOGICAL :: l_efield = .FALSE.
      LOGICAL :: l_efield_scr = .FALSE.
      LOGICAL :: l_borneffcharge = .FALSE.
      LOGICAL :: l_polar = .FALSE.
      LOGICAL :: l_symVacLevel = .TRUE. ! Symmetrize the vacua levels  

      REAL, ALLOCATABLE :: qvec(:,:)
      REAL, ALLOCATABLE :: qvec_efield(:,:)

      LOGICAL :: calcEigenVec    = .TRUE.
      
   CONTAINS
      PROCEDURE :: read_xml => read_xml_juPhon
      PROCEDURE :: mpi_bc => mpi_bc_juPhon
      PROCEDURE :: init => init_juPhon
      PROCEDURE :: precheck_juPhon
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

      !CALL mpi_bc(this%l_potout, rank, mpi_comm)
      !CALL mpi_bc(this%l_eigout, rank, mpi_comm)
      CALL mpi_bc(this%l_dfpt, rank, mpi_comm)
      CALL mpi_bc(this%l_intp, rank, mpi_comm)
      CALL mpi_bc(this%l_band, rank, mpi_comm)
      CALL mpi_bc(this%l_dos, rank, mpi_comm)
      CALL mpi_bc(this%l_scf, rank, mpi_comm)
      CALL mpi_bc(this%l_elph, rank, mpi_comm)
      CALL mpi_bc(this%l_sumrule, rank, mpi_comm)
      CALL mpi_bc(this%l_rm_qhdf, rank, mpi_comm)
      CALL mpi_bc(this%startq, rank, mpi_comm)
      CALL mpi_bc(this%stopq, rank, mpi_comm)
      CALL mpi_bc(this%qmode, rank, mpi_comm)
      CALL mpi_bc(this%i_integration, rank, mpi_comm)
      CALL mpi_bc(this%smearingGauss, rank, mpi_comm)
      CALL mpi_bc(this%qvec, rank, mpi_comm)
      CALL mpi_bc(this%qvec_efield, rank, mpi_comm)
      CALL mpi_bc(this%l_phonon, rank, mpi_comm)
      CALL mpi_bc(this%l_efield, rank, mpi_comm)
      CALL mpi_bc(this%l_efield_scr, rank, mpi_comm)
      CALL mpi_bc(this%l_borneffcharge, rank, mpi_comm)
      CALL mpi_bc(this%l_polar, rank, mpi_comm)
      CALL mpi_bc(this%qlim,rank,mpi_comm)
      CALL mpi_bc(this%gmaxzLocal,rank,mpi_comm)
      CALL mpi_bc(this%l_symVacLevel, rank, mpi_comm)
      CALL mpi_bc(this%eDiffcut, rank, mpi_comm)
      CALL mpi_bc(this%fDiffcut, rank, mpi_comm)

   END SUBROUTINE mpi_bc_juPhon

   SUBROUTINE read_xml_juPhon(this, xml)
      USE m_types_xml
      USE m_judft

      IMPLICIT NONE

      CLASS(t_juPhon), INTENT(INOUT) :: this
      TYPE(t_xml), INTENT(INOUT)     :: xml

      INTEGER::numberNodes
      CHARACTER(len=100) :: xPathA,valueString

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

        !  numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_jpTest')

        !  IF (numberNodes == 1) THEN
        !    this%l_jpTest  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_jpTest'))
        !  END IF

        !  numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_potout')

        !  IF (numberNodes == 1) THEN
        !    this%l_potout  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_potout'))
        !  END IF

        !  numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_eigout')

        !  IF (numberNodes == 1) THEN
        !    this%l_eigout  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_eigout'))
        !  END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@eDiffcut')

         IF (numberNodes == 1) THEN
           this%eDiffcut  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@eDiffcut'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@fDiffcut')

         IF (numberNodes == 1) THEN
           this%fDiffcut  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@fDiffcut'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@qlim')

         IF (numberNodes == 1) THEN
           this%qlim  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@qlim'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@gmaxzLocal')

         IF (numberNodes == 1) THEN
           this%gmaxzLocal  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@gmaxzLocal'))
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

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_elph')

         IF (numberNodes == 1) THEN
           this%l_elph  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_elph'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_sumrule')

         IF (numberNodes == 1) THEN
           this%l_sumrule  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_sumrule'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_rm_qhdf')

         IF (numberNodes == 1) THEN
           this%l_rm_qhdf  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_rm_qhdf'))
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

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@i_integration')

         IF (numberNodes == 1) THEN
           this%i_integration  = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@i_integration'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@smearingGauss')

         IF (numberNodes == 1) THEN
           this%smearingGauss  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@smearingGauss'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_phonon')

         IF (numberNodes == 1) THEN
           this%l_phonon    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_phonon'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_efield')

         IF (numberNodes == 1) THEN
           this%l_efield    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_efield'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_efield_scr')

         IF (numberNodes == 1) THEN
           this%l_efield_scr    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_efield_scr'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_borneffcharge')

         IF (numberNodes == 1) THEN
           this%l_borneffcharge    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_borneffcharge'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_polar')

         IF (numberNodes == 1) THEN
           this%l_polar    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_polar'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon/@l_symVacLevel')

         IF (numberNodes == 1) THEN
          this%l_symVacLevel    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_symVacLevel'))
         END IF


         allocate(this%qvec(0,0))
         this%qvec=xml%read_q_list('/fleurInput/output/juPhon/qVectors')
      ENDIF

      ! Before we exit check needed parameters 
      IF (this%l_dfpt) CALL this%precheck_juPhon(xml)
      
   END SUBROUTINE read_xml_juPhon

   SUBROUTINE init_juPhon(this,cell,input)

    USE m_types_cell
    USE m_types_input
    USE m_inv3
    USE m_constants

      CLASS(t_juPhon),     INTENT(INOUT) :: this
      TYPE(t_cell),    INTENT(IN)     :: cell
      TYPE(t_input),   intent(in)     :: input
      INTEGER                            :: iDir
      REAL                               :: qvec_ext(3), qvec_int(3)


      integer :: iq 
      real :: tmp_vec(3)

      if (this%l_efield .or. this%l_borneffcharge) then  
        allocate(this%qvec_efield(3,3))
        do iDir = 1,3
          qvec_ext(:) = 0.0
          qvec_int(:) = 0.0
          qvec_ext(iDir) = this%qlim
          qvec_int = matmul(qvec_ext,cell%amat)/(tpi_const) ! tpi_const is in principle irrelevant, but included for consistency with the previous q vectors
          this%qvec_efield(:,iDir) = qvec_int
        end do
      end if 

      if (input%film .and. allocated(this%qvec)) then
         if (size(this%qvec,2) > 1) THEN
            ! Due to stability we do not calculate the Gamma-Point in the case of 
            ! Film-DFPT but slighlty next to it  
            do iq = 1 , size(this%qvec,2)
               if (norm2(this%qvec(:,iq)) .LT. 1e-8 ) then 
                  tmp_vec=0.0
                  if (iq == 1 ) then ! starting q-Point --> interpolate to iQ+1
                     tmp_vec = this%qvec(:,iq+1) - this%qvec(:,iq) 
                  else ! interpolate to iQ-1 
                     tmp_vec = this%qvec(:,iq-1) - this%qvec(:,iq) 
                  end if 
                  ! construct a vector with length 1 
                  tmp_vec = tmp_vec / norm2(tmp_vec)
                  ! add to the gamma point
                  this%qvec(:,iq) = this%qvec(:,iq) + this%qlim*tmp_vec  
               end if  
            end do
         end if
      end if 


   END SUBROUTINE init_juPhon

   SUBROUTINE precheck_juPhon(this,xml)
    USE m_types_xml
    USE m_judft 

    IMPLICIT NONE 

    CLASS(t_juPhon), INTENT(IN) :: this
    TYPE(t_xml), INTENT(INOUT)  :: xml 

    INTEGER :: numberNodes
    CHARACTER(len=100) :: xPathA,valueString
    LOGICAL :: l_flag

    xPathA = '/fleurInput/calculationSetup/cutoffs/@numbands'
    numberNodes = xml%GetNumberOfNodes(xPathA)
    IF(numberNodes.EQ.1) THEN
      valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))))) 
      IF(.NOT. TRIM(ADJUSTL(valueString)).EQ.'all') CALL juDFT_error("numbands is not set to all", calledby="types_juPhon.F90")
    END IF 

   END SUBROUTINE precheck_juPhon
END MODULE m_types_juPhon
