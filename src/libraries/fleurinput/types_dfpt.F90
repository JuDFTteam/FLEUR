!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_dfpt
   USE m_judft
   USE m_types_fleurinput_base
   USE m_types_kpts
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_fleurinput_base) :: t_dfpt
      LOGICAL :: l_dfpt    = .FALSE. ! Phonon calculation on/off
      REAL    :: eDiffcut  = 1e-5   ! Cutoff for energy differences
      REAL    :: fDiffcut  = 1e-7    ! Cutoff for occupation differences
      REAL    :: qlim      = 1./100     ! qlim value
      REAL    :: gmaxzLocal   = 0.0  ! Local Gmaxz cutoff for film 

      LOGICAL :: l_intp = .FALSE.    ! Interpolate the q-set onto another one
      LOGICAL :: l_band = .FALSE.    ! Interpolate the q-set to a bandstructure
      LOGICAL :: l_dos  = .FALSE.     ! Calculate the phonon density of states
      LOGICAL :: l_scf  = .TRUE.     ! Do a self-consistency run for dynmats
      LOGICAL :: l_postprocess = .FALSE. ! Postprocessing of charge density response 
      LOGICAL :: l_elph = .FALSE.    ! Calculate electron-phonon matrix elements
      LOGICAL :: l_sumrule_scf  = .FALSE. ! Apply sumrule for dynmats in scf calculation 
      LOGICAL :: l_sumrule_intp  = .FALSE. ! Apply sumrule for dynmats in interpolation calculation 
      LOGICAL :: l_rm_qhdf  = .TRUE. ! Remove q*hdf files, after convergence
      INTEGER :: startq = 1          ! Start the q-loop at a specific point
      INTEGER :: stopq  = 0          ! Stop  the q-loop at a specific point
      INTEGER :: i_integration = 1   ! choose integration scheme for ph-linewidth (experimental)
      REAL    :: smearingGauss = 1e-7 ! Gaussian smearing for binning in el-ph
      LOGICAL :: l_phonon = .FALSE.    
      LOGICAL :: l_efield = .FALSE.
      LOGICAL :: l_efield_scr = .FALSE.
      LOGICAL :: l_borneffcharge = .FALSE.
      LOGICAL :: l_polar = .FALSE.
      LOGICAL :: l_bfield = .FALSE.
      LOGICAL :: l_symVacLevel = .TRUE. ! Symmetrize the vacua levels  
      LOGICAL :: l_WSinterpol = .TRUE.

      REAL, ALLOCATABLE :: qvec_efield(:,:)

      TYPE(t_kpts) :: qvec            ! q vectors for scf 
      TYPE(t_kpts) :: qpts_interpol   ! q mesh for interpolation
      TYPE(t_kpts) :: kpts_interpol   ! k mesh for Wannier interpolation

      INTEGER, ALLOCATABLE :: bandWindow(:)  ! Window of Blochstates we want to consider

      LOGICAL :: calcEigenVec    = .TRUE.
   CONTAINS
      PROCEDURE :: read_xml => read_xml_dfpt
      PROCEDURE :: mpi_bc => mpi_bc_dfpt
      PROCEDURE :: init => init_dfpt
      PROCEDURE :: precheck_dfpt
   END TYPE t_dfpt

   PUBLIC t_dfpt

CONTAINS

   SUBROUTINE mpi_bc_dfpt(this, mpi_comm, irank)
      USE m_mpi_bc_tool

      CLASS(t_dfpt),     INTENT(INOUT) :: this
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
      CALL mpi_bc(this%l_postprocess, rank, mpi_comm)
      CALL mpi_bc(this%l_elph, rank, mpi_comm)
      CALL mpi_bc(this%l_sumrule_scf, rank, mpi_comm)
      CALL mpi_bc(this%l_sumrule_intp, rank, mpi_comm)
      CALL mpi_bc(this%l_rm_qhdf, rank, mpi_comm)
      CALL mpi_bc(this%startq, rank, mpi_comm)
      CALL mpi_bc(this%stopq, rank, mpi_comm)
      CALL mpi_bc(this%i_integration, rank, mpi_comm)
      CALL mpi_bc(this%smearingGauss, rank, mpi_comm)
      CALL this%qvec%mpi_bc(mpi_comm, rank)
      CALL this%qpts_interpol%mpi_bc(mpi_comm, rank)
      CALL mpi_bc(this%qvec_efield, rank, mpi_comm)
      CALL mpi_bc(this%l_phonon, rank, mpi_comm)
      CALL mpi_bc(this%l_efield, rank, mpi_comm)
      CALL mpi_bc(this%l_efield_scr, rank, mpi_comm)
      CALL mpi_bc(this%l_borneffcharge, rank, mpi_comm)
      CALL mpi_bc(this%l_polar, rank, mpi_comm)
      CALL mpi_bc(this%l_bfield, rank, mpi_comm)
      CALL mpi_bc(this%qlim,rank,mpi_comm)
      CALL mpi_bc(this%gmaxzLocal,rank,mpi_comm)
      CALL mpi_bc(this%l_symVacLevel, rank, mpi_comm)
      CALL mpi_bc(this%eDiffcut, rank, mpi_comm)
      CALL mpi_bc(this%fDiffcut, rank, mpi_comm)
      CALL mpi_bc(this%bandWindow, rank, mpi_comm)
      CALL mpi_bc(this%l_WSinterpol, rank, mpi_comm)
      CALL this%kpts_interpol%mpi_bc(mpi_comm, rank)


   END SUBROUTINE mpi_bc_dfpt

   SUBROUTINE read_kpts_list(xml, listName, kpts_out)
      !! Read the k-point list <listName> from inp.xml/kpts.xml into a t_kpts
      !! (populates nkpt, bk and wtkpt).
      USE m_types_xml
      TYPE(t_xml),      INTENT(INOUT) :: xml
      CHARACTER(len=*), INTENT(IN)    :: listName
      TYPE(t_kpts),     INTENT(INOUT) :: kpts_out

      kpts_out%nkpt = xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/kPointLists/kPointList[@name="'//TRIM(listName)//'"]/kPoint')
      IF (kpts_out%nkpt > 0) THEN
         IF (ALLOCATED(kpts_out%bk))    DEALLOCATE(kpts_out%bk)
         IF (ALLOCATED(kpts_out%wtkpt)) DEALLOCATE(kpts_out%wtkpt)
         ALLOCATE(kpts_out%bk(3, kpts_out%nkpt))
         ALLOCATE(kpts_out%wtkpt(kpts_out%nkpt))
         IF (.NOT. kpts_out%read_kpts_by_name(TRIM(xml%filename_add_xml)//"inp.xml", TRIM(listName))) &
            CALL juDFT_error("Could not read kPointList '"//TRIM(listName)//"'", calledby="types_dfpt.F90")
      END IF
   END SUBROUTINE read_kpts_list

   SUBROUTINE read_qvec_list(xml, path, kpts_out)
      !! Read an inline <qVectors>/<q> list at <path> into a t_kpts
      !! (sets nkpt and fills bk; assigns uniform weights).
      USE m_types_xml
      TYPE(t_xml),      INTENT(INOUT) :: xml
      CHARACTER(len=*), INTENT(IN)    :: path
      TYPE(t_kpts),     INTENT(INOUT) :: kpts_out

      REAL, ALLOCATABLE :: q(:,:)
      INTEGER :: n

      q = xml%read_q_list(TRIM(path))     ! (3,n) from the <q> nodes
      n = SIZE(q, 2)
      IF (n > 0) THEN
         IF (ALLOCATED(kpts_out%bk))    DEALLOCATE(kpts_out%bk)
         IF (ALLOCATED(kpts_out%wtkpt)) DEALLOCATE(kpts_out%wtkpt)
         ALLOCATE(kpts_out%bk(3, n))
         ALLOCATE(kpts_out%wtkpt(n))
         kpts_out%bk    = q
         kpts_out%nkpt  = n
         kpts_out%wtkpt = 1.0 / real(n)   ! inline <q> has no weights; uniform default
      END IF
   END SUBROUTINE read_qvec_list

   SUBROUTINE read_xml_dfpt(this, xml)
      USE m_types_xml
      USE m_judft
      USE m_types_kpts

      IMPLICIT NONE

      CLASS(t_dfpt), INTENT(INOUT) :: this
      TYPE(t_xml), INTENT(INOUT)     :: xml

      INTEGER::numberNodes
      CHARACTER(len=100) :: xPathA,valueString
      CHARACTER(len=40) :: qptsListName
      TYPE(t_kpts) :: qpts_from_kpts

      REAL, ALLOCATABLE :: tmp_arr(:)


      numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt')

      IF (numberNodes == 1) THEN
         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_dfpt')

         IF (numberNodes == 1) THEN
           this%l_dfpt    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_dfpt'))
         END IF

        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_scf')

         IF (numberNodes == 1) THEN
           this%l_scf  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_scf'))
         END IF


         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_phonon')

         IF (numberNodes == 1) THEN
           this%l_phonon    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_phonon'))
         END IF

        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_efield')

         IF (numberNodes == 1) THEN
           this%l_efield    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_efield'))
         END IF

        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_borneffcharge')

         IF (numberNodes == 1) THEN
           this%l_borneffcharge    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_borneffcharge'))
         END IF

          numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_interpolate')

         IF (numberNodes == 1) THEN
           this%l_intp  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_interpolate'))
         END IF

          numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_postprocess')

         IF (numberNodes == 1) THEN
           this%l_postprocess  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_postprocess'))
         END IF

           numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_rm_qhdf')

         IF (numberNodes == 1) THEN
           this%l_rm_qhdf  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_rm_qhdf'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/@l_bfield')

         IF (numberNodes == 1) THEN
           this%l_bfield    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/@l_bfield'))
         END IF
         
        ! Phonon read-in
        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon')

        if ((numberNodes /= 1) .and. this%l_phonon) call juDFT_error("Please specify the phonon calculation with the phonon tag",calledby="types_dfpt.F90")

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/@eDiffcut')

         IF (numberNodes == 1) THEN
           this%eDiffcut  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@eDiffcut'))
         END IF 

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/@fDiffcut')

         IF (numberNodes == 1) THEN
           this%fDiffcut  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@fDiffcut'))
         END IF

          numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/@gmaxzLocal')

         IF (numberNodes == 1) THEN
           this%gmaxzLocal  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@gmaxzLocal'))
         END IF

          numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/@l_sumrule')

         IF (numberNodes == 1) THEN
           this%l_sumrule_scf  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@l_sumrule'))
         END IF

        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/@startq')

         IF (numberNodes == 1) THEN
           this%startq  = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@startq'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/@stopq')

         IF (numberNodes == 1) THEN
           this%stopq  = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@stopq'))
         END IF

          numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/phonon/@l_symVacLevel')

         IF (numberNodes == 1) THEN
          this%l_symVacLevel    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@l_symVacLevel'))
         END IF

        ! inline <qVectors> list (if present)
        call read_qvec_list(xml, '/fleurInput/output/dfpt/phonon/qVectors', this%qvec)

        ! a named kPointList overrides the inline list
         IF (xml%GetNumberOfNodes('/fleurInput/output/dfpt/phonon/@qptsListName') == 1) THEN
          qptsListName = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/dfpt/phonon/@qptsListName')))
          call read_kpts_list(xml, qptsListName, this%qvec)
         END IF


        ! efield read-in
        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/efield')

        if ( (numberNodes /= 1) .and. this%l_efield) call juDFT_error("Please specify the phonon calculation with the phonon tag",calledby="types_dfpt.F90")

        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/efield/@l_efield_scr')

         IF (numberNodes == 1) THEN
           this%l_efield_scr    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/efield/@l_efield_scr'))
         END IF

        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/efield/@qlim')

         IF (numberNodes == 1) THEN
           this%qlim  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/efield/@qlim'))
         END IF


        ! interpolation read-in 
        numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/interpolation')

        if ((numberNodes /= 1) .and. this%l_intp) call juDFT_error("Please specify the interpolation with the interpolation tag",calledby="types_dfpt.F90")
         
         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/interpolation/@l_band')

         IF (numberNodes == 1) THEN
           this%l_band  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/interpolation/@l_band'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/interpolation/@l_dos')

         IF (numberNodes == 1) THEN
           this%l_dos  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/interpolation/@l_dos'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/interpolation/@l_WSinterpol')

         IF (numberNodes == 1) THEN
          this%l_WSinterpol    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/interpolation/@l_WSinterpol'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/interpolation/@l_polar')

         IF (numberNodes == 1) THEN
           this%l_polar    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/interpolation/@l_polar'))
         END IF

          numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/interpolation/@l_sumrule')

         IF (numberNodes == 1) THEN
           this%l_sumrule_intp  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/interpolation/@l_sumrule'))
         END IF


          IF (xml%GetNumberOfNodes('/fleurInput/output/dfpt/interpolation/@qptsListName') == 1) THEN
          qptsListName = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/dfpt/interpolation/@qptsListName')))
          call read_kpts_list(xml, qptsListName, this%qpts_interpol)
        END IF

        ! postprocess read-in
         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/postprocess')
         ! introduce switch for postprocessing 
        if ( (numberNodes /= 1) .and. this%l_postprocess) call juDFT_error("Please specify the interpolation with the interpolation tag",calledby="types_dfpt.F90")
        

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/postprocess/@l_elph')

         IF (numberNodes == 1) THEN
           this%l_elph  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/postprocess/@l_elph'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/postprocess/@i_integration')

         IF (numberNodes == 1) THEN
           this%i_integration  = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/postprocess/@i_integration'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/postprocess/@smearingGauss')

         IF (numberNodes == 1) THEN
           this%smearingGauss  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/dfpt/postprocess/@smearingGauss'))
         END IF

         numberNodes = xml%GetNumberOfNodes('/fleurInput/output/dfpt/postprocess/@bandWindow')

         IF (numberNodes == 1) THEN
          allocate(this%bandWindow(2))
          allocate(tmp_arr(2))
          valueString = xml%GetAttributeValue('/fleurInput/output/dfpt/postprocess/@bandWindow')
          call evaluateList(tmp_arr,valueString)
          this%bandWindow = tmp_arr
         END IF


          IF (xml%GetNumberOfNodes('/fleurInput/output/dfpt/postprocess/@qMeshName') == 1) THEN
          qptsListName = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/dfpt/postprocess/@qMeshName')))
          call read_kpts_list(xml, qptsListName, this%qpts_interpol)
        END IF

         IF (xml%GetNumberOfNodes('/fleurInput/output/dfpt/postprocess/@kMeshName') == 1) THEN
          qptsListName = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/dfpt/postprocess/@kMeshName')))
          call read_kpts_list(xml, qptsListName, this%kpts_interpol)
        END IF


      END IF

      ! Before we exit check needed parameters 
      IF (this%l_dfpt) CALL this%precheck_dfpt(xml)
      
   END SUBROUTINE read_xml_dfpt

   SUBROUTINE init_dfpt(this,cell,input,sym,noco)

    USE m_types_cell
    USE m_types_input
    USE m_inv3
    USE m_constants
    USE m_types_sym
    USE m_types_noco

      CLASS(t_dfpt), INTENT(INOUT) :: this
      TYPE(t_cell),    INTENT(IN)     :: cell
      TYPE(t_input),   intent(in)     :: input
      TYPE(t_sym),     intent(in)     :: sym
      TYPE(t_noco),    intent(in)     :: noco
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

      ! if (input%film .and. allocated(this%qvec)) then
      !    if (size(this%qvec,2) > 1) THEN
      !       ! Due to stability we do not calculate the Gamma-Point in the case of 
      !       ! Film-DFPT but slighlty next to it  
      !       do iq = 1 , size(this%qvec,2)
      !          if (norm2(this%qvec(:,iq)) .LT. 1e-8 ) then 
      !             tmp_vec=0.0
      !             if (iq == 1 ) then ! starting q-Point --> interpolate to iQ+1
      !                tmp_vec = this%qvec(:,iq+1) - this%qvec(:,iq) 
      !             else ! interpolate to iQ-1 
      !                tmp_vec = this%qvec(:,iq-1) - this%qvec(:,iq) 
      !             end if 
      !             ! construct a vector with length 1 
      !             tmp_vec = tmp_vec / norm2(tmp_vec)
      !             ! add to the gamma point
      !             this%qvec(:,iq) = this%qvec(:,iq) + this%qlim*tmp_vec  
      !          end if  
      !       end do
      !    end if
      ! end if 

 
      ! Build the full BZ (nkptf, bkf, bkp, bksym) for the interpolation meshes.
      if (allocated(this%qvec%bk)) then 
        call this%qvec%init(sym, input%film, .false., (noco%l_soc.or.noco%l_ss))
        this%qvec%nkpt3 = this%qvec%calcNkpt3()
      end if 
      if (allocated(this%qpts_interpol%bk)) then 
        call this%qpts_interpol%init(sym, input%film, .false., (noco%l_soc.or.noco%l_ss))
        this%qpts_interpol%nkpt3 = this%qpts_interpol%calcNkpt3()
      end if 
      if (allocated(this%kpts_interpol%bk)) then 
        call this%kpts_interpol%init(sym, input%film, .false., (noco%l_soc.or.noco%l_ss))
        this%kpts_interpol%nkpt3 = this%kpts_interpol%calcNkpt3() 
      end if 
   END SUBROUTINE init_dfpt

   SUBROUTINE precheck_dfpt(this,xml)
    ! This routine checks if the general input structure 
    ! for the dfpt calculation is compatible 
    USE m_types_xml
    USE m_judft 

    IMPLICIT NONE 

    CLASS(t_dfpt), INTENT(IN) :: this
    TYPE(t_xml), INTENT(INOUT)  :: xml 

    INTEGER :: numberNodes
    CHARACTER(len=100) :: xPathA,valueString
    LOGICAL :: l_flag

    xPathA = '/fleurInput/calculationSetup/cutoffs/@numbands'
    numberNodes = xml%GetNumberOfNodes(xPathA)
    IF(numberNodes.EQ.1) THEN
      valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))))) 
      IF(.NOT. TRIM(ADJUSTL(valueString)).EQ.'all') CALL juDFT_error("numbands is not set to all", calledby="types_dfpt.F90")
    END IF 

    if (this%l_phonon ) then 
      if (allocated(this%qvec%bk)) then
        if (this%qvec%nkpt .eq. 0 ) call juDFT_warn("No q-Points were given while trying to do a phonon calculation. Please insert q points",calledby="types_dfpt.F90")
      end if
    end if 

   END SUBROUTINE precheck_dfpt
END MODULE m_types_dfpt
