!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hub1inp
   USE m_juDFT
   USE m_constants
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_fleurinput_base):: t_hub1inp
      !Convergence criteria for the density matrix
      INTEGER :: itmax = 5
      REAL    :: minoccDistance=1.0e-2
      REAL    :: minmatDistance=1.0e-3
      LOGICAL :: l_dftspinpol=.FALSE.     !Determines whether the DFT part is spin-polarized in a magnetic DFT+Hubbard 1 calculation
      LOGICAL :: l_fullMatch=.TRUE.       !Determines whether two chemical potentials are used to match (if possible)
      LOGICAL :: l_nonsphDC=.TRUE.        !Determines whether to remove the nonspherical contributions to the Hamiltonian (in the HIA orbital)
      LOGICAL :: l_correctEtot = .TRUE.   !Perform additional scf cycle without spin averaging of the correlated shell in DFT with frozen density matrix
      LOGICAL :: l_forceHIAiteration = .FALSE.

      !Parameters for the solver
      REAL     :: beta = 100.0 !inverse temperature
      INTEGER  :: n_occpm = 2  !number of particle excitations considered in the solver

      REAL, ALLOCATABLE :: init_occ(:) !initial occupation
      REAL, ALLOCATABLE :: ccf(:) !crystal field factor
      REAL, ALLOCATABLE :: cfCoeffs(:,:,:)
      REAL, ALLOCATABLE :: xi_par(:) !Fixed SOC parameters

      INTEGER, ALLOCATABLE :: n_exc(:)
      INTEGER, ALLOCATABLE :: exc_l(:,:) !l quantum number from which the intraorbital exchange
      REAL,    ALLOCATABLE :: exc(:,:) !exchange splitting parameter
      REAL,    ALLOCATABLE :: init_mom(:,:) !initial magnetic moment

      !Additional arguments to be passed on to hloc.cfg (at the moment only real)
      INTEGER,             ALLOCATABLE :: n_addArgs(:)
      CHARACTER(len=100),  ALLOCATABLE :: arg_keys(:,:)
      REAL,                ALLOCATABLE :: arg_vals(:,:)

      CHARACTER(len=100), allocatable :: post_process_tasks(:)

      !Switches for arguments that were explicitly given and should not be calculated from DFT
      LOGICAL,ALLOCATABLE :: l_soc_given(:)
      LOGICAL,ALLOCATABLE :: l_ccf_given(:)

   CONTAINS
      PROCEDURE :: read_xml   => read_xml_hub1inp
      PROCEDURE :: mpi_bc     => mpi_bc_hub1inp
   END TYPE t_hub1inp
   PUBLIC t_hub1inp

CONTAINS

   SUBROUTINE mpi_bc_hub1inp(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_hub1inp), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF
      CALL mpi_bc(this%itmax,rank,mpi_comm)
      CALL mpi_bc(this%minoccDistance,rank,mpi_comm)
      CALL mpi_bc(this%minmatDistance,rank,mpi_comm)
      CALL mpi_bc(this%l_dftspinpol,rank,mpi_comm)
      CALL mpi_bc(this%l_fullMatch,rank,mpi_comm)
      CALL mpi_bc(this%l_nonsphDC,rank,mpi_comm)
      CALL mpi_bc(this%l_forceHIAiteration,rank,mpi_comm)
      CALL mpi_bc(this%beta,rank,mpi_comm)
      CALL mpi_bc(this%n_occpm,rank,mpi_comm)
      CALL mpi_bc(this%init_occ,rank,mpi_comm)
      CALL mpi_bc(this%ccf,rank,mpi_comm)
      CALL mpi_bc(this%cfCoeffs,rank,mpi_comm)
      CALL mpi_bc(this%xi_par,rank,mpi_comm)
      CALL mpi_bc(this%n_exc,rank,mpi_comm)
      CALL mpi_bc(this%exc_l,rank,mpi_comm)
      CALL mpi_bc(this%exc,rank,mpi_comm)
      CALL mpi_bc(this%init_mom,rank,mpi_comm)
      CALL mpi_bc(this%n_addArgs,rank,mpi_comm)
      CALL mpi_bc(this%arg_keys,rank,mpi_comm)
      CALL mpi_bc(this%arg_vals,rank,mpi_comm)
      CALL mpi_bc(this%l_soc_given,rank,mpi_comm)
      CALL mpi_bc(this%l_ccf_given,rank,mpi_comm)
      CALL mpi_bc(this%post_process_tasks,rank,mpi_comm)
   END SUBROUTINE mpi_bc_hub1inp

   SUBROUTINE read_xml_hub1inp(this, xml)
      USE m_types_xml
      CLASS(t_hub1inp), INTENT(INOUT):: this
      TYPE(t_xml),INTENT(INOUT) ::xml

      INTEGER::numberNodes,ntype,n_maxaddArgs, numTasks
      INTEGER::i_hia,itype,i_exc,i_addArg,i,j,hub1_l,i_cf,l,m
      CHARACTER(len=100)  :: xPathA,xPathB,xPathS,key,tmp_str
      CHARACTER(len=300)  :: tasks
      REAL::val

      ntype = xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      n_maxaddArgs = 5 !Maximum allowed number of additional arguments (excluding xiSOC and ccf)

      ALLOCATE(this%init_occ(4*ntype),source=0.0)
      ALLOCATE(this%ccf(4*ntype),source=1.0)
      ALLOCATE(this%cfCoeffs(4*ntype,0:6,-6:6),source=0.0)
      ALLOCATE(this%xi_par(4*ntype),source=0.001)
      ALLOCATE(this%n_exc(4*ntype),source=0)
      ALLOCATE(this%exc_l(4*ntype,lmaxU_const),source=-1)
      ALLOCATE(this%exc(4*ntype,lmaxU_const),source=0.0)
      ALLOCATE(this%init_mom(4*ntype,lmaxU_const),source=0.0)
      ALLOCATE(this%n_addArgs(4*ntype),source=0)
      ALLOCATE(this%arg_keys(4*ntype,n_maxaddArgs))
      ALLOCATE(this%post_process_tasks(4))
      this%arg_keys='' !For some reason source doesn't work here
      this%post_process_tasks='' !For some reason source doesn't work here
      ALLOCATE(this%arg_vals(4*ntype,n_maxaddArgs),source=0.0)
      ALLOCATE(this%l_soc_given(4*ntype),source=.FALSE.)
      ALLOCATE(this%l_ccf_given(4*ntype),source=.FALSE.)

      !General parameters:
      xPathA = '/fleurInput/calculationSetup/ldaHIA'
      numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
      IF(numberNodes==1) THEN
         this%itmax = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@itmaxHubbard1'))
         this%minoccDistance = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minoccDistance'))
         this%minmatDistance = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minmatDistance'))
         this%beta = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@beta'))
         this%n_occpm = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n_occpm'))
         this%l_dftspinpol = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dftspinpol'))
         this%l_fullMatch = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@fullMatch'))
         this%l_nonsphDC = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_nonsphDC'))
         this%l_correctEtot = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_correctEtot'))
         IF(xml%versionNumber>=34) THEN
            this%l_forceHIAiteration = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_forceHIAiteration'))
         ENDIF

         if(xml%versionNumber>=36) then
            if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/postProcess') == 1) then
               tasks = xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/postProcess/text()')
               numTasks = xml%countStringTokens(tasks)
               deallocate(this%post_process_tasks)
               allocate(this%post_process_tasks(numTasks))
               do i = 1, numTasks
                  this%post_process_tasks(i) = xml%popFirstStringToken(tasks)
               enddo
            endif
         endif
      ENDIF

      !Read in the additional information given in the ldaHIA tags (exchange splitting and additional keywords)
      i_hia=0
      DO itype = 1, ntype
         xPathS = xml%speciesPath(itype)
         DO j = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/ldaHIA')
            i_hia = i_hia + 1
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/ldaHIA[',j,']'
            !Read in the hubbard 1 orbital for a later check
            hub1_l = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l'))

            !Initial occupation
            tmp_str = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@init_occ')))
            IF(TRIM(ADJUSTL(tmp_str))=="calc") THEN
               this%init_occ(i_hia) = -9e99
            ELSE
               this%init_occ(i_hia) = evaluateFirstOnly(TRIM(ADJUSTL(tmp_str)))
            ENDIF

            !Additional exchange splitting
            DO i_exc = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/exc')
               WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'/exc[',i_exc,']'
               IF(i_exc>lmaxU_const) CALL juDFT_error("Too many additional exchange splittings provided. Maximum is 3.",&
                                           calledby="read_xml_hub1inp")
               this%n_exc(i_hia) = this%n_exc(i_hia) + 1
               this%exc_l(i_hia,i_exc) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l'))
               this%exc(i_hia,i_exc) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@J'))
               tmp_str = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@init_mom')))
               IF(TRIM(ADJUSTL(tmp_str))=="calc") THEN
                  this%init_mom(i_hia,i_exc) = -9e99
               ELSE
                  this%init_mom(i_hia,i_exc) = evaluateFirstOnly(TRIM(ADJUSTL(tmp_str)))
               ENDIF

               !Check if the given l is valid (l<3 and not the same as the hubbard orbital)
               IF(this%exc_l(i_hia,i_exc).EQ.hub1_l.OR.this%exc_l(i_hia,i_exc).GT.3) &
                  CALL juDFT_error("Additional exchange splitting: Not a valid l"&
                                  ,calledby="read_xml_hub1inp")
               !Check if there already is a defined exchange splitting on this orbital
               DO i = 1, this%n_exc(i_hia)-1
                  IF(this%exc_l(i_hia,i_exc)==this%exc_l(i_hia,i)) &
                     CALL juDFT_error("Two exchange splittings defined for equal l"&
                                     ,calledby="read_xml_hub1inp")
               ENDDO
            ENDDO

            DO i_cf = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/cFCoeff')
               WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'/cFCoeff[',i_cf,']'

               l = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l'))
               m = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@m'))
               !Check if the given l is valid (l<3 and not the same as the hubbard orbital)
               IF(l.LT.0 .OR. l.GT.6 .OR. ABS(m).GT.l) &
                  CALL juDFT_error("Crystal Field Coefficient: Not a valid l,m combination (|m|<l and l<=6)"&
                                  ,calledby="read_xml_hub1inp")

               IF(ABS(this%cfCoeffs(i_hia,l,m)).GT.1e-12) &
                  CALL juDFT_error("Crystal Field Coefficient: Two Coefficients for the same lm"&
                                  ,calledby="read_xml_hub1inp")
               this%cfCoeffs(i_hia,l,m) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@value'))
               this%cfCoeffs(i_hia,l,-m) = this%cfCoeffs(i_hia,l,m)
            ENDDO

            DO i_addArg = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/addArg')

               WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'/addArg[',i_addArg,']'
               IF(i_addArg>n_maxaddArgs) CALL juDFT_error("Too many additional arguments provided. Maximum is 5.",&
                                                          calledby="read_xml_hub1inp")

               key = xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@key')
               val = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@value'))

               DO i = 1, this%n_addArgs(i_hia)
                  IF(TRIM(ADJUSTL(key)).EQ.TRIM(ADJUSTL(this%arg_keys(i_hia,i)))) THEN
                     CALL juDFT_error("Ambigous additional arguments: You specified two arguments with the same keyword"&
                                     ,calledby="read_xml_hub1inp")
                  ENDIF
               ENDDO

               SELECT CASE(TRIM(ADJUSTL(key)))
               CASE('xiSOC')
                  !Do not get soc from DFT and use provided value
                  IF(this%l_soc_given(i_hia)) CALL juDFT_error("Two SOC parameters provided",calledby="read_xml_hub1inp")
                  this%l_soc_given(i_hia) = .TRUE.
                  this%xi_par(i_hia) = val
                  IF(ABS(this%xi_par(i_hia))< 0.001)  this%xi_par(i_hia) = 0.001
               CASE('ccf')
                  IF(this%l_ccf_given(i_hia)) CALL juDFT_error("Two crystal field factors provided",calledby="read_xml_hub1inp")
                  this%l_ccf_given(i_hia) = .TRUE.
                  this%ccf(i_hia) = val
               CASE DEFAULT
                  !Additional argument -> simply pass on to solver
                  this%n_addArgs(i_hia) = this%n_addArgs(i_hia) + 1
                  this%arg_keys(i_hia,this%n_addArgs(i_hia)) = TRIM(ADJUSTL(key))
                  this%arg_vals(i_hia,this%n_addArgs(i_hia)) = val
               END SELECT
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE read_xml_hub1inp

END MODULE m_types_hub1inp
