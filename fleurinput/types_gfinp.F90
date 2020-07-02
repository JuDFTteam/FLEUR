!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_gfinp
   USE m_juDFT
   USE m_types_fleurinput_base
   USE m_constants

   IMPLICIT NONE
   PRIVATE

   !Define constants for the shape argument
   INTEGER, PARAMETER :: CONTOUR_RECTANGLE_CONST  = 1
   INTEGER, PARAMETER :: CONTOUR_SEMICIRCLE_CONST = 2
   INTEGER, PARAMETER :: CONTOUR_DOS_CONST        = 3

   TYPE t_gfelementtype
      SEQUENCE
      !defines the l and atomType elements for given greens function element
      !(used for mapping index in types_greensf)
      INTEGER :: l = -1
      INTEGER :: lp = -1
      INTEGER :: atomType = 0
      INTEGER :: atomTypep = 0
      INTEGER :: iContour = 0 !Which energy contour is used
      REAL    :: atomDiff(3) = [0.0,0.0,0.0] !Distance between atoms for intersite phase
      LOGICAL :: l_sphavg = .TRUE. !Is this element calculated with or without radial dependence
      !Parameters for the determination of the upper cutoff of the Kramers-Kronig Integration
      LOGICAL :: l_fixedCutoffset = .FALSE.
      REAL    :: fixedCutoff = 0.0
      INTEGER :: refCutoff = -1 !Choose cutoff to be the same as another GFelement
   END TYPE t_gfelementtype

   TYPE t_contourInp
      !Contains the input parameters for a contour
      INTEGER :: shape = 2 !If no contour is specified write out the standard Semicircle contour (for inpgen)
      CHARACTER(len=100) :: label = "default" !name of the contour
      !Endpoints of the complex energy contour
      REAL    :: eb = -1.0
      REAL    :: et = 0.0
      !Shape-specific parameters (for details look at e_contour in types_greensf.f90)
      !Parameters for mode == 1 (Rectangle)
      INTEGER :: n1 = 10
      INTEGER :: n2 = 128
      INTEGER :: n3 = 20
      INTEGER :: nmatsub = 5
      REAL    :: sigma = 0.005
      !Parameters for mode == 2 (Semicircle)
      INTEGER :: ncirc = 128
      REAL    :: alpha = 1.0
      !Switches for mode == 3 (Equidistant mesh shifted by sigmaDOS)
      INTEGER :: nDOS = 1301
      REAL    :: sigmaDOS = 0.0314
      LOGICAL :: l_anacont = .FALSE. !Determines wether to include an analytical continuation at the edges
      LOGICAL :: l_dosfermi = .FALSE.!Determines wether the integration weights include the fermi distribution
   END TYPE t_contourInp

   TYPE, EXTENDS(t_fleurinput_base):: t_gfinp
      !General logical switches
      LOGICAL :: l_mperp = .FALSE.
      REAL    :: minCalcDistance=-1.0 !This distance has to be reached before green's functions are calculated
                                      !Negative means it is evaluated at every iteration
      !Number of elements
      INTEGER :: n = 0
      !Information on the elements to be calculated
      TYPE(t_gfelementtype), ALLOCATABLE :: elem(:)
      !Parameters for the energy mesh on the real axis
      INTEGER :: ne    = 2700
      REAL    :: ellow = -1.0
      REAL    :: elup  =  1.0
      INTEGER :: numberContours = 0
      TYPE(t_contourInp), ALLOCATABLE :: contour(:)

      !Arrays to indicate that certain Green's Functions are used for special calculations
      INTEGER, ALLOCATABLE :: hiaElem(:)
      INTEGER, ALLOCATABLE :: torgueElem(:,:)
      INTEGER, ALLOCATABLE :: numTorgueElems(:)
   CONTAINS
      PROCEDURE :: read_xml       => read_xml_gfinp
      PROCEDURE :: mpi_bc         => mpi_bc_gfinp
      PROCEDURE :: init           => init_gfinp
      PROCEDURE :: find           => find_gfelem
      PROCEDURE :: find_contour   => find_contour
      PROCEDURE :: add            => add_gfelem
      PROCEDURE :: uniqueElements => uniqueElements_gfinp
      PROCEDURE :: eMesh          => eMesh_gfinp
      PROCEDURE :: checkRadial    => checkRadial_gfinp
      PROCEDURE :: checkSphavg    => checkSphavg_gfinp
      PROCEDURE :: addNearestNeighbours => addNearestNeighbours_gfelem
   END TYPE t_gfinp

   PUBLIC t_gfinp, t_contourInp, t_gfelementtype
   PUBLIC CONTOUR_RECTANGLE_CONST, CONTOUR_SEMICIRCLE_CONST, CONTOUR_DOS_CONST

CONTAINS

   SUBROUTINE mpi_bc_gfinp(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_gfinp), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank,myrank,n,ierr
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF
      CALL mpi_bc(this%l_mperp,rank,mpi_comm)
      CALL mpi_bc(this%minCalcDistance,rank,mpi_comm)
      CALL mpi_bc(this%n,rank,mpi_comm)
      CALL mpi_bc(this%ne,rank,mpi_comm)
      CALL mpi_bc(this%ellow,rank,mpi_comm)
      CALL mpi_bc(this%elup,rank,mpi_comm)
      CALL mpi_bc(this%numberContours,rank,mpi_comm)
      CALL mpi_bc(this%hiaElem,rank,mpi_comm)
      CALL mpi_bc(this%torgueElem,rank,mpi_comm)
      CALL mpi_bc(this%numTorgueElems,rank,mpi_comm)

#ifdef CPP_MPI
      CALL mpi_COMM_RANK(mpi_comm,myrank,ierr)
      IF (myrank.NE.rank) THEN
         IF (ALLOCATED(this%elem)) DEALLOCATE(this%elem)
         IF (ALLOCATED(this%contour)) DEALLOCATE(this%contour)
         ALLOCATE(this%elem(this%n))
         ALLOCATE(this%contour(this%numberContours))
      ENDIF
      DO n=1,this%n
         CALL mpi_bc(this%elem(n)%l,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%atomType,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%lp,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%atomTypep,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%iContour,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%l_fixedCutoffset,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%fixedCutoff,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%refCutoff,rank,mpi_comm)
         CALL mpi_bc(rank,mpi_comm,this%elem(n)%atomDiff)
         CALL mpi_bc(this%elem(n)%l_sphavg,rank,mpi_comm)
      ENDDO
      DO n=1,this%numberContours
         CALL mpi_bc(this%contour(n)%shape,rank,mpi_comm)
         CALL mpi_bc(rank,mpi_comm,this%contour(n)%label) !Order reversed for fixed routine
         CALL mpi_bc(this%contour(n)%eb,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%et,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%n1,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%n2,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%n3,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%nmatsub,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%sigma,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%ncirc,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%alpha,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%nDOS,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%sigmaDOS,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%l_anacont,rank,mpi_comm)
         CALL mpi_bc(this%contour(n)%l_dosfermi,rank,mpi_comm)
      ENDDO
#endif

   END SUBROUTINE mpi_bc_gfinp

   SUBROUTINE read_xml_gfinp(this, xml)
      USE m_types_xml
      CLASS(t_gfinp), INTENT(INOUT):: this
      TYPE(t_xml),INTENT(INOUT) ::xml

      INTEGER :: numberNodes,ntype,itype,n_hia,i_gf,refL,refGF,nshells
      INTEGER :: i,l,lp,iContour,iContourp
      REAL    :: fixedCutoff
      CHARACTER(len=200)  :: xPathA,xPathS,label,cutoffArg,str
      CHARACTER(len=1),PARAMETER :: spdf(0:3) = ['s','p','d','f']
      LOGICAL :: l_gfinfo_given,l_fixedCutoffset,l_sphavg
      LOGICAL :: lp_calc(0:3,0:3)

      xPathA = '/fleurInput/calculationSetup/greensFunction'
      numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
      l_gfinfo_given = numberNodes.EQ.1

      IF (l_gfinfo_given) THEN
         this%l_mperp=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
         this%minCalcDistance=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minCalcDistance'))

         xPathA = '/fleurInput/calculationSetup/greensFunction/realAxis'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF(numberNodes /= 1) CALL juDFT_error("Error reading in gf-information: realAxis not specified correctly",&
                                               calledby="read_xml_gfinp")

         this%ne = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ne'))
         this%ellow = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ellow'))
         this%elup = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@elup'))

         !Read in the complex energy contours

         !How many defined contours are there
         xPathA = '/fleurInput/calculationSetup/greensFunction/contourRectangle'
         this%numberContours = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         xPathA = '/fleurInput/calculationSetup/greensFunction/contourSemicircle'
         this%numberContours = this%numberContours + xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         xPathA = '/fleurInput/calculationSetup/greensFunction/contourDOS'
         this%numberContours = this%numberContours + xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))

         IF(this%numberContours.EQ.0) THEN
            CALL juDFT_error("Error reading in gf-information: No complex energy contour specified",&
                             calledby="read_xml_gfinp")
         ENDIF

         ALLOCATE(this%contour(this%numberContours))

         iContour = 0
         xPathS = '/fleurInput/calculationSetup/greensFunction'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/contourRectangle')
         DO i = 1, numberNodes
            iContour = iContour + 1
            WRITE(xPathA,'(a,i0,a)') TRIM(ADJUSTL(xPathS))//'/contourRectangle[',i,']'
            this%contour(iContour)%shape = CONTOUR_RECTANGLE_CONST
            this%contour(iContour)%eb = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eb'))
            !et cannot be varied from the fermi energy for this contour
            this%contour(iContour)%n1 = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n1'))
            this%contour(iContour)%n2 = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n2'))
            this%contour(iContour)%n3 = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n3'))
            this%contour(iContour)%nmatsub = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nmatsub'))
            this%contour(iContour)%sigma = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sigma'))
            this%contour(iContour)%label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
         ENDDO

         xPathS = '/fleurInput/calculationSetup/greensFunction'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/contourSemicircle')
         DO i = 1, numberNodes
            iContour = iContour + 1
            WRITE(xPathA,'(a,i0,a)') TRIM(ADJUSTL(xPathS))//'/contourSemicircle[',i,']'
            this%contour(iContour)%shape = CONTOUR_SEMICIRCLE_CONST
            this%contour(iContour)%eb = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eb'))
            this%contour(iContour)%et = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@et'))
            this%contour(iContour)%ncirc = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n'))
            this%contour(iContour)%alpha = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@alpha'))
            this%contour(iContour)%label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
         ENDDO

         xPathS = '/fleurInput/calculationSetup/greensFunction'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/contourDOS')
         DO i = 1, numberNodes
            iContour = iContour + 1
            WRITE(xPathA,'(a,i0,a)') TRIM(ADJUSTL(xPathS))//'/contourDOS[',i,']'
            this%contour(iContour)%shape = CONTOUR_DOS_CONST
            this%contour(iContour)%eb = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eb'))
            this%contour(iContour)%et = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@et'))
            this%contour(iContour)%nDOS = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n'))
            this%contour(iContour)%sigmaDOS = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sigma'))
            this%contour(iContour)%l_anacont = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@analytical_cont'))
            this%contour(iContour)%l_dosfermi = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_fermi'))
            this%contour(iContour)%label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
         ENDDO
         !Check for ambigous labels
         DO iContour = 1, this%numberContours
            DO iContourp = iContour+1, this%numberContours
               IF(TRIM(ADJUSTL(this%contour(iContour)%label)) == TRIM(ADJUSTL(this%contour(iContourp)%label))) THEN
                  CALL juDFT_error("Ambigous definition of energy contours", calledby="read_xml_gfinp")
               ENDIF
            ENDDO
         ENDDO

      ENDIF

      ntype = xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      n_hia = 0

      ALLOCATE(this%elem((lmaxU_const+1)**2*ntype))
      ALLOCATE(this%hiaElem(4*ntype))
      ALLOCATE(this%numTorgueElems(ntype),source=0)
      ALLOCATE(this%torgueElem(ntype,(lmaxU_const+1)**2),source=-1)

      DO itype = 1, ntype
         xPathS=xml%speciesPath(itype)

         !Read in all possible tags, which need a greens function calculation

         !Declaration of a general Green's Function Calculation
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/greensfCalculation')
            WRITE(xPathA,'(a,i0,a)') TRIM(ADJUSTL(xPathS))//'/greensfCalculation[',i,']'

            label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
            l_sphavg = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_sphavg'))
            iContour = this%find_contour(TRIM(ADJUSTL(label)))
            cutoffArg = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kkintgrCutoff')))
            nshells = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nshells'))

            refL = -1
            SELECT CASE(TRIM(ADJUSTL(cutoffArg)))
            CASE('calc')
               !calculate the cutoff from the number of states
               l_fixedCutoffset = .FALSE.
            CASE('s','p','d','f')
               !Reference cutoff given (will be processed after)
               l_fixedCutoffset = .FALSE.
               IF(TRIM(ADJUSTL(cutoffArg))=='s') refL = 0
               IF(TRIM(ADJUSTL(cutoffArg))=='p') refL = 1
               IF(TRIM(ADJUSTL(cutoffArg))=='d') refL = 2
               IF(TRIM(ADJUSTL(cutoffArg))=='f') refL = 3
            CASE default
               !Fixed cutoff set
               fixedCutoff = evaluateFirstOnly(TRIM(ADJUSTL(cutoffArg)))
               l_fixedCutoffset = .TRUE.
            END SELECT

            numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/matrixElements')
            IF(numberNodes==1) THEN
               xPathA = TRIM(ADJUSTL(xPathA))//'/matrixElements'
               DO l = 0,lmaxU_const
                  str = xml%GetAttributeValue(TRIM(xPathA)//'/@'//spdf(l))
                  READ(str,'(4l2)') (lp_calc(lp,l),lp=0,3)
                  DO lp = 0,lmaxU_const
                     IF(.NOT.lp_calc(lp,l)) CYCLE
                     i_gf =  this%add(l,itype,iContour,l_sphavg,lp=lp,l_fixedCutoffset=l_fixedCutoffset,&
                                   fixedCutoff=fixedCutoff,l_inter=(nshells/=0))
                  ENDDO
               ENDDO
            ENDIF

            numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/diagElements')
            IF(numberNodes==1) THEN
               lp_calc = .FALSE.
               xPathA = TRIM(ADJUSTL(xPathA))//'/diagElements'
               DO l = 0,lmaxU_const
                  lp_calc(l,l) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@'//spdf(l)))
                  IF(.NOT.lp_calc(l,l)) CYCLE
                  i_gf =  this%add(l,itype,iContour,l_sphavg,l_fixedCutoffset=l_fixedCutoffset,&
                                   fixedCutoff=fixedCutoff,l_inter=(nshells/=0))
               ENDDO
            ENDIF

            !Find the reference element
            IF(refL /= -1) THEN
               !Find the element
               refGF = this%find(refL,itype,iContour,l_sphavg)
               DO l = 0,lmaxU_const
                  DO lp = 0,lmaxU_const
                     IF(.NOT.lp_calc(lp,l)) CYCLE
                     i_gf = this%find(l,itype,iContour,l_sphavg,lp=lp)
                     IF(i_gf==refGF) CYCLE
                     this%elem(i_gf)%refCutoff = refGF
                  ENDDO
               ENDDO
            ENDIF

         ENDDO
         !Declaration of a DFT+Hubbard 1 calculation
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/ldaHIA')
            WRITE(xPathA,'(a,i0,a)') TRIM(ADJUSTL(xPathS))//'/ldaHIA[',i,']'
            !No offdiagonal l-part (only a single element)
            l = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l'))
            label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
            iContour = this%find_contour(TRIM(ADJUSTL(label)))
            cutoffArg = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kkintgrCutoff')))
            IF(TRIM(ADJUSTL(cutoffArg))=="calc") THEN
               l_fixedCutoffset = .FALSE.
            ELSE
               fixedCutoff = evaluateFirstOnly(TRIM(ADJUSTL(cutoffArg)))
               l_fixedCutoffset = .TRUE.
            ENDIF
            !Hubbard 1 GF has to be spherically averaged
            i_gf =  this%add(l,itype,iContour,.TRUE.,l_fixedCutoffset=l_fixedCutoffset,&
                             fixedCutoff=fixedCutoff)
            n_hia = n_hia + 1
            this%hiaElem(n_hia) = i_gf
         ENDDO

         WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/torgueCalculation'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF(numberNodes==1) THEN
            label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
            iContour = this%find_contour(TRIM(ADJUSTL(label)))
            cutoffArg = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kkintgrCutoff')))

            refL = -1
            SELECT CASE(TRIM(ADJUSTL(cutoffArg)))
            CASE('calc')
               !calculate the cutoff from the number of states
               l_fixedCutoffset = .FALSE.
            CASE('s','p','d','f')
               !Reference cutoff given (will be processed after)
               l_fixedCutoffset = .FALSE.
               IF(TRIM(ADJUSTL(cutoffArg))=='s') refL = 0
               IF(TRIM(ADJUSTL(cutoffArg))=='p') refL = 1
               IF(TRIM(ADJUSTL(cutoffArg))=='d') refL = 2
               IF(TRIM(ADJUSTL(cutoffArg))=='f') refL = 3
            CASE default
               !Fixed cutoff set
               fixedCutoff = evaluateFirstOnly(TRIM(ADJUSTL(cutoffArg)))
               l_fixedCutoffset = .TRUE.
            END SELECT

            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/torgueCalculation/greensfElements'
            DO l = 0,lmaxU_const
               str = xml%GetAttributeValue(TRIM(xPathA)//'/'//spdf(l))
               READ(str,'(4l2)') (lp_calc(lp,l),lp=0,3)
               DO lp = 0,lmaxU_const
                  IF(.NOT.lp_calc(lp,l)) CYCLE
                  !Torgue GF has to have radial dependence
                  i_gf =  this%add(l,itype,iContour,.FALSE.,lp=lp,l_fixedCutoffset=l_fixedCutoffset,&
                                   fixedCutoff=fixedCutoff)
                  this%numTorgueElems(itype) = this%numTorgueElems(itype) + 1
                  this%torgueElem(itype,this%numTorgueElems(itype)) = i_gf
               ENDDO
            ENDDO

            !Find the reference element
            IF(refL /= -1) THEN
               !Find the element
               refGF = this%find(refL,itype,iContour,.FALSE.)
               DO l = 0,lmaxU_const
                  DO lp = 0,lmaxU_const
                     IF(.NOT.lp_calc(lp,l)) CYCLE
                     i_gf =  this%find(l,itype,iContour,.FALSE.,lp=lp)
                     IF(i_gf==refGF) CYCLE
                     this%elem(i_gf)%refCutoff = refGF
                  ENDDO
               ENDDO
            ENDIF

         ENDIF

      ENDDO

      IF(this%n>0 .AND. .NOT.l_gfinfo_given) THEN
         CALL juDFT_error("Error reading in gf-information: No general information found for the gf-calculations",&
                           calledby="read_xml_gfinp")
      ENDIF

      !Check the input for validity
      IF(this%n.GT.0) THEN
         IF(this%elup.GT.1.0) CALL juDFT_warn("Cutoff for the Greens function calculation should never be higher"//&
                                              "than 1htr above efermi",calledby="read_xml_gfinp")
         IF(this%elup.LE.this%ellow) CALL juDFT_error("Not a valid energy grid elup<ellow",calledby="read_xml_gfinp")
         IF(ANY(this%elem(:this%n)%l.LT.2)) CALL juDFT_warn("Green's function for s and p orbitals not tested",&
                                                            calledby="read_xml_gfinp")
         IF(ANY(this%elem(:this%n)%l.GT.3)) CALL juDFT_error("Green's function only implemented for l<=3",&
                                                             calledby="read_xml_gfinp")
      ENDIF

   END SUBROUTINE read_xml_gfinp

   SUBROUTINE init_gfinp(this,atoms,sym,noco,input,l_write)

      USE m_types_atoms
      USE m_types_sym
      USE m_types_noco
      USE m_types_input

      CLASS(t_gfinp),   INTENT(INOUT)  :: this
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_input),    INTENT(IN)     :: input
      LOGICAL,          INTENT(IN)     :: l_write

      INTEGER :: i_gf,l,lp,atomType,atomTypep,iContour
      LOGICAL :: l_inter,l_offd,l_sphavg,l_interAvg,l_offdAvg
      INTEGER :: hiaElem(atoms%n_hia)

      IF(this%n==0) RETURN !Nothing to do here

      !Find the elements for which we need to compute the nearest neighbours
      DO i_gf = 1, this%n
         l  = this%elem(i_gf)%l
         lp = this%elem(i_gf)%lp
         atomType  = this%elem(i_gf)%atomType
         atomTypep = this%elem(i_gf)%atomTypep
         iContour = this%elem(i_gf)%iContour

         IF(atomTypep==-1) THEN
            !Replace the current element by the onsite one
            this%elem(i_gf)%atomTypep = atomType
            CALL this%addNearestNeighbours(1,l,lp,iContour,atomType,this%elem(i_gf)%l_fixedCutoffset,&
                                           this%elem(i_gf)%fixedCutoff,this%elem(i_gf)%refCutoff,atoms,sym)
         ENDIF
      ENDDO

      !Reallocate with correct size
      hiaElem = this%hiaElem(:atoms%n_hia)
      IF(ALLOCATED(this%hiaElem)) DEALLOCATE(this%hiaElem)
      ALLOCATE(this%hiaElem(atoms%n_hia))
      this%hiaElem = hiaElem

      !Input checks
      IF(l_write) THEN
         IF(this%l_mperp.AND..NOT.noco%l_mperp) THEN
            CALL juDFT_error("For l_mperp for Green's Functions the l_mperp switch for noco has to be True",&
                             calledby="init_gfinp")
         ENDIF

         l_inter = .FALSE.
         l_offd = .FALSE.
         l_interAvg = .FALSE.
         l_offdAvg = .FALSE.
         DO i_gf = 1, this%n
            l  = this%elem(i_gf)%l
            lp = this%elem(i_gf)%lp
            atomType  = this%elem(i_gf)%atomType
            atomTypep = this%elem(i_gf)%atomTypep
            l_sphavg  = this%elem(i_gf)%l_sphavg
            IF(atomType.NE.atomTypep) THEN
               l_inter = .TRUE.
               IF(l_sphavg) l_interAvg = .TRUE.
            ENDIF
            IF(l.NE.lp) THEN
               l_offd = .TRUE.
               IF(l_sphavg) l_offdAvg = .TRUE.
            ENDIF

         ENDDO

         IF(l_inter) THEN
            IF(sym%nop>1) THEN
                  CALL juDFT_error("Symmetries and intersite Green's Function not implemented",&
                                   calledby="init_gfinp")
            ELSE IF(l_interAvg) THEN
               CALL juDFT_error("Spherical average and intersite Green's Function not implemented",&
                                calledby="init_gfinp")
            ENDIF
         ENDIF

         IF(l_offd) THEN
            IF(sym%nop>1) THEN
               CALL juDFT_warn("Symmetries and l-offdiagonal Green's Function not correctly implemented",&
                                calledby="init_gfinp")
            ELSE IF(l_offdAvg) THEN
               CALL juDFT_error("Spherical average and l-offdiagonal Green's Function not implemented",&
                                calledby="init_gfinp")
            ENDIF
         ENDIF

         IF(this%minCalcDistance>=0.0) THEN
            IF(input%mindistance>this%minCalcDistance) THEN
               CALL juDFT_warn("The minimum Distance for Green's Function Calculation"// &
                               "is smaller than the distance requirement:"//&
                               "No Green's Functions will be calculated", calledby="init_gfinp")
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE init_gfinp

   FUNCTION uniqueElements_gfinp(this,ind,l_sphavg,indUnique) Result(uniqueElements)

      CLASS(t_gfinp),   INTENT(IN)     :: this
      INTEGER, OPTIONAL,INTENT(IN)     :: ind
      LOGICAL, OPTIONAL,INTENT(IN)     :: l_sphavg !uniqueElements are determined separately for radial dependence and spherically averaging
      INTEGER, OPTIONAL,INTENT(INOUT)  :: indUnique !Position of the corresponding unique Element for a given ind

      INTEGER :: uniqueElements !Number of unique elements before ind or in the whole array (if ind is not present)

      INTEGER :: maxGF
      INTEGER :: l,lp,atomType,atomTypep,iUnique,iContour,i_gf
      LOGICAL :: l_sphavgArg, l_sphavgElem

      l_sphavgArg = .TRUE.
      IF(PRESENT(l_sphavg)) l_sphavgArg = l_sphavg

      uniqueElements = 0

      IF(PRESENT(ind)) THEN
         maxGF = ind
      ELSE
         maxGF = this%n
      ENDIF
      DO i_gf = 1, maxGF
         l  = this%elem(i_gf)%l
         lp = this%elem(i_gf)%lp
         atomType  = this%elem(i_gf)%atomType
         atomTypep = this%elem(i_gf)%atomTypep
         iContour  = this%elem(i_gf)%iContour
         l_sphavgElem  = this%elem(i_gf)%l_sphavg
         IF(l_sphavgElem .neqv. l_sphavgArg) CYCLE
         iUnique   = this%find(l,atomType,iContour,l_sphavgElem,lp=lp,nTypep=atomTypep,&
                               uniqueMax=i_gf)

         IF(iUnique == i_gf) uniqueElements = uniqueElements +1
      ENDDO

      IF(PRESENT(indUnique)) THEN
         IF(.NOT.PRESENT(ind)) CALL juDFT_error("ind and indUnique have to be provided at the same time",&
                                                calledby="uniqueElements_gfinp")
         l  = this%elem(ind)%l
         lp = this%elem(ind)%lp
         atomType  = this%elem(ind)%atomType
         atomTypep = this%elem(ind)%atomTypep
         iContour  = this%elem(ind)%iContour
         l_sphavgElem = this%elem(ind)%l_sphavg

         indUnique = this%find(l,atomType,iContour,l_sphavgElem,lp=lp,nTypep=atomTypep,&
                               uniqueMax=ind)
      ENDIF

   END FUNCTION uniqueElements_gfinp

   INTEGER FUNCTION add_gfelem(this,l,nType,iContour,l_sphavg,lp,nTypep,l_fixedCutoffset,fixedCutoff,l_inter) Result(i_gf)

      CLASS(t_gfinp),      INTENT(INOUT)  :: this
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: nType
      INTEGER,             INTENT(IN)     :: iContour
      LOGICAL,             INTENT(IN)     :: l_sphavg
      INTEGER, OPTIONAL,   INTENT(IN)     :: lp
      INTEGER, OPTIONAL,   INTENT(IN)     :: nTypep !Specify the second atom
      LOGICAL, OPTIONAL,   INTENT(IN)     :: l_fixedCutoffset
      REAL,    OPTIONAL,   INTENT(IN)     :: fixedCutoff
      LOGICAL, OPTIONAL,   INTENT(IN)     :: l_inter!To be used in init when atoms is not available and nTypep was not specified


      LOGICAL l_found

      IF(PRESENT(l_inter).AND.PRESENT(nTypep)) CALL juDFT_error("Conflicting arguments: l_inter and nTypep given",&
                                                                hint="This is a bug in FLEUR, please report",&
                                                                calledby="add_gfelem")

      !Check if this job has already been added
      i_gf = this%find(l,nType,iContour,l_sphavg,lp=lp,nTypep=nTypep,l_found=l_found)
      IF(l_found) RETURN !Element was found

      this%n = this%n + 1
      i_gf = this%n
      this%elem(this%n)%l = l
      this%elem(this%n)%atomType = nType
      this%elem(this%n)%iContour = iContour
      this%elem(this%n)%l_sphavg = l_sphavg
      IF(PRESENT(lp)) THEN
         this%elem(this%n)%lp = lp
      ELSE
         this%elem(this%n)%lp = l
      ENDIF
      IF(PRESENT(l_inter)) THEN
         IF(l_inter) THEN
            !Temporary index to mark later in gfinp%init
            this%elem(this%n)%atomTypep = -1
         ELSE
            this%elem(this%n)%atomTypep = nType
         ENDIF
      ELSE IF(PRESENT(nTypep)) THEN
         !Explicit declaration of element
         this%elem(this%n)%atomTypep = nTypep
      ELSE
         !No intersite element
         this%elem(this%n)%atomTypep = nType
      ENDIF
      IF(PRESENT(l_fixedCutoffset)) THEN
         IF(.NOT.PRESENT(fixedCutoff)) CALL juDFT_error("l_fixedCutoffset Present without fixedCutoff", &
                                                        hint="This is a bug in FLEUR please report",&
                                                        calledby="add_gfelem")
         this%elem(this%n)%l_fixedCutoffset = l_fixedCutoffset
         IF(l_fixedCutoffset) THEN
            this%elem(this%n)%fixedCutoff = fixedCutoff
         ENDIF
      ENDIF

   END FUNCTION add_gfelem

   SUBROUTINE addNearestNeighbours_gfelem(this,nshells,l,lp,refAtom,iContour,l_fixedCutoffset,fixedCutoff,refCutoff,atoms,sym)

      USE m_types_atoms
      USE m_types_sym

      !TODO: atoms outside the unit cell

      CLASS(t_gfinp),   INTENT(INOUT)  :: this
      INTEGER,          INTENT(IN)     :: nshells !How many nearest neighbour shells are requested
      INTEGER,          INTENT(IN)     :: l
      INTEGER,          INTENT(IN)     :: lp
      INTEGER,          INTENT(IN)     :: refAtom !which is the reference atom
      INTEGER,          INTENT(IN)     :: iContour
      LOGICAL,          INTENT(IN)     :: l_fixedCutoffset
      REAL,             INTENT(IN)     :: fixedCutoff
      INTEGER,          INTENT(IN)     :: refCutoff
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym

      INTEGER :: ishell,natomp,i_gf
      REAL :: minDist
      REAL, ALLOCATABLE :: dist(:)

      IF(sym%nop>1) CALL juDFT_error("nearest neighbour GF + symmetries not implemented",&
                                     calledby="addNearestNeighbours_gfelem")

      ALLOCATE(dist(atoms%nat),source=0.0)
      DO natomp = 1, atoms%nat
         dist(natomp) = SQRT((atoms%taual(1,natomp) - atoms%taual(1,refAtom))**2 + &
                             (atoms%taual(2,natomp) - atoms%taual(2,refAtom))**2 + &
                             (atoms%taual(3,natomp) - atoms%taual(3,refAtom))**2)
         IF(ABS(dist(natomp)).LT.1e-12) dist(natomp) = 9e99 !The atom itself was already added
      ENDDO

      ishell = 0
      DO WHILE(ishell<=nshells)
         !search for the atoms with the current minimal distance
         minDist = MINVAL(dist)
         DO natomp = 1, atoms%nat
            IF(ABS(dist(natomp)-minDist).LT.1e-12) THEN
               !Add the element to the gfinp%elem array
               this%elem(this%n)%atomDiff = atoms%taual(:,natomp) - atoms%taual(:,refAtom)
               !cannot be spherically averaged
               i_gf =  this%add(l,refAtom,iContour,.FALSE.,lp=lp,nTypep=natomp,l_fixedCutoffset=l_fixedCutoffset,&
                                fixedCutoff=fixedCutoff)
               this%elem(i_gf)%refCutoff = refCutoff
               dist(natomp) = 9e99 !Eliminate from the list
            ENDIF
         ENDDO
         ishell = ishell + 1
      ENDDO

   END SUBROUTINE addNearestNeighbours_gfelem

   INTEGER FUNCTION find_gfelem(this,l,nType,iContour,l_sphavg,lp,nTypep,uniqueMax,l_found) result(i_gf)

      !Maps between the four indices (l,lp,nType,nTypep) and the position in the
      !gf arrays

      CLASS(t_gfinp),      INTENT(IN)    :: this
      INTEGER,             INTENT(IN)    :: l
      INTEGER,             INTENT(IN)    :: nType
      INTEGER,             INTENT(IN)    :: iContour
      LOGICAL,             INTENT(IN)    :: l_sphavg
      INTEGER, OPTIONAL,   INTENT(IN)    :: lp
      INTEGER, OPTIONAL,   INTENT(IN)    :: nTypep

      INTEGER, OPTIONAL,   INTENT(IN)    :: uniqueMax  !If uniqueMax is present it will return the
                                                       !index of the unique element, meaning
                                                       !the same (l,lp,type,typep) but different contours

      LOGICAL, OPTIONAL,   INTENT(INOUT) :: l_found    !If this switch is not provided the program
                                                       !will assume that the element has to be present and
                                                       !terminate with an error message if the
                                                       !element is not found (for adding elements)

      LOGICAL :: search

      search = .TRUE.
      i_gf = 0

      DO WHILE(search)
         i_gf = i_gf + 1

         IF(i_gf>this%n) THEN
            !Did not find the element
            IF(PRESENT(l_found)) THEN
               i_gf = -1
               search = .FALSE.
               l_found = .FALSE.
               CYCLE
            ELSE
               CALL juDFT_error("Green's function element not found",&
                                hint="This is a bug in FLEUR, please report",&
                                calledby="find_gfelem")
            ENDIF
         ENDIF
         !--------------------------------------------
         ! Check the current element
         !--------------------------------------------
         IF(this%elem(i_gf)%l.NE.l) CYCLE
         IF(this%elem(i_gf)%atomType.NE.nType) CYCLE
         IF(PRESENT(lp)) THEN
            IF(this%elem(i_gf)%lp.NE.lp) CYCLE
         ELSE
            IF(this%elem(i_gf)%lp.NE.l) CYCLE
         ENDIF
         IF(PRESENT(nTypep)) THEN
            IF(nTypep/=nType) THEN
               IF(this%elem(i_gf)%atomTypep.NE.nTypep) CYCLE
            ELSE
               !The -1 will be replaced with the onsite element in init_gfinp
               IF(this%elem(i_gf)%atomTypep.NE.nTypep.AND.this%elem(i_gf)%atomTypep.NE.-1) CYCLE
            ENDIF
         ELSE
            !The -1 will be replaced with the onsite element in init_gfinp
            IF(this%elem(i_gf)%atomTypep.NE.nType.AND.this%elem(i_gf)%atomTypep.NE.-1) CYCLE
         ENDIF
         IF(this%elem(i_gf)%l_sphavg .neqv. l_sphavg) CYCLE
         !If we are here and smaller than uniqueMax the element is not unique
         !i.e they only differ in the choice of the energy contour
         IF(PRESENT(uniqueMax)) THEN
            IF(i_gf>uniqueMax) CALL juDFT_error('i_gf>uniqueMax',calledby="find_gfelem")
            RETURN
         ENDIF
         IF(this%elem(i_gf)%iContour.NE.iContour) CYCLE

         !If we are here we found the element
         IF(PRESENT(l_found)) l_found=.TRUE.
         search = .FALSE.
      ENDDO
   END FUNCTION find_gfelem

   FUNCTION find_contour(this,label) result(iContour)

      !Finds the contour defined with the given label

      CLASS(t_gfinp),      INTENT(IN)  :: this
      CHARACTER(len=*),  INTENT(IN)  :: label

      INTEGER :: iContour
      LOGICAL :: search

      search = .TRUE.
      iContour = 0
      DO WHILE(search)
         iContour = iContour + 1

         IF(iContour>this%numberContours) THEN
               CALL juDFT_error("Energy contour not found",&
                                hint="Check the labels of your contours",&
                                calledby="find_contour")
         ENDIF
         !--------------------------------------------
         ! Check the current element
         !--------------------------------------------
         IF(TRIM(ADJUSTL(this%contour(iContour)%label)).NE.TRIM(ADJUSTL(label))) CYCLE
         search = .FALSE.
      ENDDO

   END FUNCTION find_contour

   SUBROUTINE eMesh_gfinp(this,ef,del,eb,eMesh)

      !Gives back the information for the energy mesh on the real axis
      !Energies are shifted according to the fermi level

      CLASS(t_gfinp),               INTENT(IN)    :: this
      REAL,                         INTENT(IN)    :: ef
      REAL, OPTIONAL,               INTENT(INOUT) :: del
      REAL, OPTIONAL,               INTENT(INOUT) :: eb
      REAL, ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: eMesh(:)

      INTEGER :: ie
      REAL :: delTmp,ebTmp

      ebTmp  = ef+this%ellow
      delTmp = (this%elup-this%ellow)/REAL(this%ne-1)

      IF(PRESENT(eb)) eb = ebTmp
      IF(PRESENT(del)) del = delTmp

      IF(PRESENT(eMesh)) THEN
         IF(ALLOCATED(eMesh)) DEALLOCATE(eMesh)
         ALLOCATE(eMesh(this%ne))

         DO ie = 1, this%ne
            eMesh(ie) = (ie-1)*delTmp+ebTmp
         ENDDO
      ENDIF

   END SUBROUTINE eMesh_gfinp

   PURE LOGICAL FUNCTION checkRadial_gfinp(this)

      !Check if there are any elements with radial dependence
      CLASS(t_gfinp),               INTENT(IN)    :: this

      INTEGER :: i_gf

      checkRadial_gfinp = .FALSE.
      DO i_gf = 1, this%n
         IF(.NOT.this%elem(i_gf)%l_sphavg) checkRadial_gfinp = .TRUE.
      ENDDO

   END FUNCTION checkRadial_gfinp

   PURE LOGICAL FUNCTION checkSphavg_gfinp(this)

      !Check if there are any elements with  spherical averaging
      CLASS(t_gfinp),               INTENT(IN)    :: this

      INTEGER :: i_gf

      checkSphavg_gfinp = .FALSE.
      DO i_gf = 1, this%n
         IF(this%elem(i_gf)%l_sphavg) checkSphavg_gfinp = .TRUE.
      ENDDO

   END FUNCTION checkSphavg_gfinp

END MODULE m_types_gfinp
