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
      !defines the l and atomType elements for given greens function element
      !(used for mapping index in gfinp%elem)
      INTEGER :: l  = -1
      INTEGER :: lp = -1
      INTEGER :: atomType  = 0
      INTEGER :: atomTypep = 0
      REAL    :: atomDiff(3)  = [0.0,0.0,0.0] !Distance between atoms (lattice coordinates) for intersite phase

      INTEGER :: iContour = 0      !Which energy contour is used
      LOGICAL :: l_sphavg = .TRUE. !Is this element calculated with or without radial dependence

      !Parameters for the determination of the upper cutoff of the Kramers-Kronig Integration
      LOGICAL :: l_fixedCutoffset = .FALSE.
      REAL    :: fixedCutoff = 0.0
      INTEGER :: refCutoff   = -1 !Choose cutoff to be the same as another
   CONTAINS
      PROCEDURE :: countLOs   => countLOs_gfelem !Count the local orbitals attached to the element
      PROCEDURE :: isoffDiag  => isOffDiag_gfelem !Is this element offdiagonal (i.e either l/=lp or intersite)
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
      LOGICAL :: l_mperp         = .FALSE.
      LOGICAL :: l_resolvent     = .FALSE.
      LOGICAL :: l_outputSphavg  = .FALSE.
      LOGICAL :: l_intFullRadial = .FALSE.
      REAL    :: minCalcDistance = -1.0 !This distance has to be reached before green's functions are calculated
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
      PROCEDURE :: read_xml             => read_xml_gfinp
      PROCEDURE :: mpi_bc               => mpi_bc_gfinp
      PROCEDURE :: init                 => init_gfinp
      PROCEDURE :: find                 => find_gfelem
      PROCEDURE :: find_contour         => find_contour
      PROCEDURE :: add                  => add_gfelem
      PROCEDURE :: addNearestNeighbours => addNearestNeighbours_gfelem
      PROCEDURE :: uniqueElements       => uniqueElements_gfinp
      PROCEDURE :: eMesh                => eMesh_gfinp

      !Checks to see if specific elements are present in the elem array
      PROCEDURE :: checkRadial          => checkRadial_gfinp !With Radial dependence
      PROCEDURE :: checkSphavg          => checkSphavg_gfinp !Without Radial dependence
      PROCEDURE :: checkOnsite          => checkOnsite_gfinp !Onsite Element (atomDiff=0 l=lp atom=atomp)
      PROCEDURE :: checkOffdiagonal     => checkOffdiagonal_gfinp !Offdiagonal Element (not Onsite)
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
      CALL mpi_bc(this%l_resolvent,rank,mpi_comm)
      CALL mpi_bc(this%minCalcDistance,rank,mpi_comm)
      CALL mpi_bc(this%l_outputSphavg,rank,mpi_comm)
      CALL mpi_bc(this%l_intFullRadial,rank,mpi_comm)
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

      INTEGER :: numberNodes,ntype,itype,n_hia,i_gf,refL,refGF,nshells,maxIntersiteElem
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
         this%l_resolvent=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_resolvent'))
         this%minCalcDistance=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minCalcDistance'))
         this%l_outputSphavg=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@outputSphavg'))
         this%l_intFullRadial=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@intFullRadial'))


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

      !Maximum number of intersite elements (number of unit cells in a 7x7x7 block) (will be reallocated later with actual size)
      maxIntersiteElem = 343
      ALLOCATE(this%elem((lmaxU_const+1)**2*maxIntersiteElem*ntype)) !ntype shoul be nat ??
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
                  str = xml%GetAttributeValue(TRIM(xPathA)//'/'//spdf(l))
                  READ(str,'(4l2)') (lp_calc(lp,l),lp=0,3)
                  DO lp = 0,lmaxU_const
                     IF(.NOT.lp_calc(lp,l)) CYCLE
                     i_gf =  this%add(l,itype,iContour,l_sphavg,lp=lp,l_fixedCutoffset=l_fixedCutoffset,&
                                   fixedCutoff=fixedCutoff,nshells=nshells)
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
                                   fixedCutoff=fixedCutoff,nshells=nshells)
               ENDDO
            ENDIF

            !Find the reference element
            IF(refL /= -1 .AND.ANY(lp_calc)) THEN
               !Find the element
               refGF = this%find(refL,itype,iContour,l_sphavg,nTypep=-nshells)
               DO l = 0,lmaxU_const
                  DO lp = 0,lmaxU_const
                     IF(.NOT.lp_calc(lp,l)) CYCLE
                     i_gf = this%find(l,itype,iContour,l_sphavg,lp=lp,nTypep=-nshells)
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
            IF(refL /= -1 .AND.ANY(lp_calc)) THEN
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

   SUBROUTINE init_gfinp(this,atoms,sym,noco,cell,input)

      USE m_types_atoms
      USE m_types_sym
      USE m_types_noco
      USE m_types_input
      USE m_types_cell

      CLASS(t_gfinp),   INTENT(INOUT)  :: this
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_cell),     INTENT(IN)     :: cell
      TYPE(t_input),    INTENT(IN)     :: input

      INTEGER :: i_gf,l,lp,atomType,atomTypep,iContour,refCutoff
      LOGICAL :: l_inter,l_offd,l_sphavg,l_interAvg,l_offdAvg
      INTEGER :: hiaElem(atoms%n_hia)
      LOGICAL :: written(atoms%nType)
      REAL    :: atomDiff(3)
      TYPE(t_gfelementtype), ALLOCATABLE :: gfelem(:)

      IF(this%n==0) RETURN !Nothing to do here

      written = .FALSE.
      !Find the elements for which we need to compute the nearest neighbours
      DO i_gf = 1, this%n
         l  = this%elem(i_gf)%l
         lp = this%elem(i_gf)%lp
         atomType  = this%elem(i_gf)%atomType
         atomTypep = this%elem(i_gf)%atomTypep
         iContour = this%elem(i_gf)%iContour
         refCutoff = this%elem(i_gf)%refCutoff
         refCutoff = MERGE(i_gf,refCutoff,refCutoff==-1) !If no refCutoff is set for the intersite element
                                                         !we take the onsite element as reference

         IF(atomTypep<0) THEN !This indicates that the nshells argument was written here
            !Replace the current element by the onsite one
            this%elem(i_gf)%atomTypep = atomType
            CALL this%addNearestNeighbours(ABS(atomTypep),l,lp,atomType,iContour,this%elem(i_gf)%l_fixedCutoffset,&
                                           this%elem(i_gf)%fixedCutoff,refCutoff,&
                                           atoms,cell,sym,.NOT.written(atomType))
            written(atomType) = .TRUE.
         ENDIF
      ENDDO

      !After this point there are no new green's function elements to be added

      !Reallocate with correct size
      hiaElem = this%hiaElem(:atoms%n_hia)
      IF(ALLOCATED(this%hiaElem)) DEALLOCATE(this%hiaElem)
      ALLOCATE(this%hiaElem(atoms%n_hia))
      this%hiaElem = hiaElem

      ALLOCATE(gfelem(this%n))
      gfelem = this%elem(:this%n)
      IF(ALLOCATED(this%elem)) DEALLOCATE(this%elem)
      ALLOCATE(this%elem(this%n))
      this%elem = gfelem

      !Input checks
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
         atomDiff  = this%elem(i_gf)%atomDiff
         IF(atomType.NE.atomTypep.OR.ANY(ABS(atomDiff).GT.1e-12)) THEN
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
               CALL juDFT_warn("Symmetries and intersite Green's Function not correctly implemented",&
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

      IF(ANY(this%numTorgueElems(:)>0)) THEN
         IF(input%jspins.NE.2) CALL juDFT_error("Torgue calculation only for magnetic systems", calledby="init_gfinp")
         IF(sym%nop>1) CALL juDFT_warn("Torgue calculation only without symmetries", calledby="init_gfinp")
      ENDIF

#ifdef CPP_DEBUG
      WRITE(*,*) "Green's Function Elements: "
      WRITE(*,'(8(A,tr5))') "l","lp","atomType","atomTypep","iContour","l_sphavg","refCutoff","atomDiff"
      DO i_gf = 1, this%n
         WRITE(*,'(5I10,1l5,I10,3f14.8)') this%elem(i_gf)%l,this%elem(i_gf)%lp,this%elem(i_gf)%atomType,this%elem(i_gf)%atomTypep,&
                                          this%elem(i_gf)%iContour,this%elem(i_gf)%l_sphavg,this%elem(i_gf)%refCutoff,&
                                          this%elem(i_gf)%atomDiff(:)
      ENDDO
#endif

   END SUBROUTINE init_gfinp

   FUNCTION uniqueElements_gfinp(this,atoms,ind,l_sphavg,lo,indUnique,maxLO) Result(uniqueElements)

      USE m_types_atoms

      CLASS(t_gfinp),   INTENT(IN)     :: this
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      INTEGER, OPTIONAL,INTENT(IN)     :: ind
      LOGICAL, OPTIONAL,INTENT(IN)     :: l_sphavg !uniqueElements are determined separately for radial dependence and spherically averaging
      LOGICAL, OPTIONAL,INTENT(IN)     :: lo       !We are interested in unique LO elems (radial dependence)
      INTEGER, OPTIONAL,INTENT(INOUT)  :: indUnique !Position of the corresponding unique Element for a given ind
      INTEGER, OPTIONAL,INTENT(INOUT)  :: maxLO    !Maximum number of Elements associated with a GF element

      INTEGER :: uniqueElements !Number of unique elements before ind or in the whole array (if ind is not present)

      INTEGER :: maxGF,nLO
      INTEGER :: l,lp,atomType,atomTypep,iUnique,iContour,i_gf
      LOGICAL :: l_sphavgArg, l_sphavgElem,loArg
      REAL    :: atomDiff(3)

      !Process optional switches and arguments
      l_sphavgArg = .TRUE.
      IF(PRESENT(l_sphavg)) l_sphavgArg = l_sphavg

      loArg = .FALSE.
      IF(PRESENT(lo)) loArg = lo

      uniqueElements = 0
      IF(PRESENT(maxLO)) maxLO = 0

      IF(PRESENT(ind)) THEN
         maxGF = ind
      ELSE
         maxGF = this%n
      ENDIF

      !Count the unique Elements before maxGF
      DO i_gf = 1, maxGF
         l  = this%elem(i_gf)%l
         lp = this%elem(i_gf)%lp
         atomType  = this%elem(i_gf)%atomType
         atomTypep = this%elem(i_gf)%atomTypep
         iContour  = this%elem(i_gf)%iContour
         l_sphavgElem  = this%elem(i_gf)%l_sphavg
         atomDiff(:) = this%elem(i_gf)%atomDiff(:)

         IF(l_sphavgElem .neqv. l_sphavgArg) CYCLE
         iUnique   = this%find(l,atomType,iContour,l_sphavgElem,lp=lp,nTypep=atomTypep,&
                               atomDiff=atomDiff,uniqueMax=i_gf)

         IF(iUnique == i_gf) THEN
            IF(loArg) THEN
               IF(.NOT.l_sphavgElem) THEN
                  nLO = this%elem(i_gf)%countLOs(atoms)
                  IF(nLO/=0) uniqueElements = uniqueElements +1
                  IF(PRESENT(maxLO)) THEN
                     IF(nLO>maxLO) maxLO = nLO
                  ENDIF
               ENDIF
            ELSE
               uniqueElements = uniqueElements +1
            ENDIF
         ENDIF
      ENDDO

      !Get the index of the unque element associated with the ind-th element
      IF(PRESENT(indUnique)) THEN
         IF(.NOT.PRESENT(ind)) CALL juDFT_error("ind and indUnique have to be provided at the same time",&
                                                calledby="uniqueElements_gfinp")
         l  = this%elem(ind)%l
         lp = this%elem(ind)%lp
         atomType  = this%elem(ind)%atomType
         atomTypep = this%elem(ind)%atomTypep
         iContour  = this%elem(ind)%iContour
         l_sphavgElem = this%elem(ind)%l_sphavg
         atomDiff(:) = this%elem(ind)%atomDiff(:)

         indUnique = this%find(l,atomType,iContour,l_sphavgElem,lp=lp,nTypep=atomTypep,&
                               atomDiff=atomDiff,uniqueMax=ind)
      ENDIF

   END FUNCTION uniqueElements_gfinp

   INTEGER FUNCTION add_gfelem(this,l,nType,iContour,l_sphavg,lp,nTypep,atomDiff,l_fixedCutoffset,fixedCutoff,nshells) Result(i_gf)

      CLASS(t_gfinp),      INTENT(INOUT)  :: this
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: nType
      INTEGER,             INTENT(IN)     :: iContour
      LOGICAL,             INTENT(IN)     :: l_sphavg
      INTEGER, OPTIONAL,   INTENT(IN)     :: lp
      INTEGER, OPTIONAL,   INTENT(IN)     :: nTypep !Specify the second atom
      REAL,    OPTIONAL,   INTENT(IN)     :: atomDiff(:)
      LOGICAL, OPTIONAL,   INTENT(IN)     :: l_fixedCutoffset
      REAL,    OPTIONAL,   INTENT(IN)     :: fixedCutoff
      INTEGER, OPTIONAL,   INTENT(IN)     :: nshells


      LOGICAL l_found

      IF(PRESENT(nshells).AND.PRESENT(nTypep)) CALL juDFT_error("Conflicting arguments: nshells and nTypep given",&
                                                                hint="This is a bug in FLEUR, please report",&
                                                                calledby="add_gfelem")

      l_found = .FALSE.
      !Check if this job has already been added
      IF(PRESENT(nTypep).OR..NOT.PRESENT(nshells)) THEN
         i_gf = this%find(l,nType,iContour,l_sphavg,lp=lp,nTypep=nTypep,atomDiff=atomDiff,l_found=l_found)
      ELSE IF(PRESENT(nshells)) THEN
         !Make sure we have not added the reference element for intersite elements
         i_gf = this%find(l,nType,iContour,l_sphavg,lp=lp,nTypep=-nshells,atomDiff=atomDiff,l_found=l_found)
      ENDIF
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
      IF(PRESENT(nshells)) THEN
         IF(nshells/=0) THEN
            !Temporary index to mark later in gfinp%init
            this%elem(this%n)%atomTypep = -nshells
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
      IF(PRESENT(atomDiff)) THEN
         this%elem(this%n)%atomDiff(:) = atomDiff(:)
      ELSE
         this%elem(this%n)%atomDiff(:) = 0.0
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

   SUBROUTINE addNearestNeighbours_gfelem(this,nshells,l,lp,refAtom,iContour,l_fixedCutoffset,fixedCutoff,&
                                          refCutoff,atoms,cell,sym,l_write)

      USE m_types_atoms
      USE m_types_cell
      USE m_types_sym
      USE m_sort
      USE m_inv3

      !This is essentially a simplified version of chkmt, because we have a given
      !reference atom and do not need to consider all distances between all atoms

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
      TYPE(t_cell),     INTENT(IN)     :: cell
      TYPE(t_sym),      INTENT(IN)     :: sym
      LOGICAL,          INTENT(IN)     :: l_write

      REAL,    PARAMETER :: tol = 1e-7

      INTEGER :: i,j,k,m,n,na,iAtom,maxAtoms,identicalAtoms,nshellDist,cubeStartIndex,cubeEndIndex
      INTEGER :: numNearestNeighbors,ishell,lastIndex,iNeighborAtom,i_gf
      INTEGER :: iop,ishell1,ishellAtom,nshellAtom,nshellAtom1,nshellsFound
      REAL :: currentDist,minDist,amatAuxDet,lastDist
      REAL :: amatAux(3,3), invAmatAux(3,3)
      REAL :: taualAux(3,atoms%nat), posAux(3,atoms%nat)
      REAL :: refPos(3),point(3),pos(3),diff(3)
      REAL :: currentDiff(3),offsetPos(3),diffRot(3)

      INTEGER, ALLOCATABLE :: nearestNeighbors(:)
      INTEGER, ALLOCATABLE :: neighborAtoms(:)
      INTEGER, ALLOCATABLE :: distIndexList(:)
      REAL,    ALLOCATABLE :: nearestNeighborDists(:)
      REAL,    ALLOCATABLE :: nearestNeighborDiffs(:,:)
      REAL,    ALLOCATABLE :: neighborAtomsDiff(:,:)
      REAL,    ALLOCATABLE :: sqrDistances(:)

      REAL,    ALLOCATABLE :: shellDistance(:)
      REAL,    ALLOCATABLE :: shellDiff(:,:,:)
      INTEGER, ALLOCATABLE :: shellAtom(:)
      INTEGER, ALLOCATABLE :: numshellAtoms(:)
      REAL,    ALLOCATABLE :: shellAux(:,:)
      REAL,    ALLOCATABLE :: shellAux1(:,:)

      CALL timestart("Green's Function: Add nearest Neighbors")
!     1. For the 1st version the auxiliary unit cell is just a copy of the original unit cell with
!        all atoms within the cell.

      !How many unit cells are included (determine based on nshells ??)
      cubeEndIndex   =  1
      cubeStartIndex = -cubeEndIndex
      !How many atoms are there in total (in all of the copies)
      maxAtoms = atoms%nat * (cubeEndIndex-cubeStartIndex+1)**3

      DO i = 1, 3
         DO j = 1, 3
            amatAux(i,j) = cell%amat(i,j)
         END DO
      END DO

      DO i = 1, atoms%nat
         taualAux(1,i) = atoms%taual(1,i) - FLOOR(atoms%taual(1,i))
         taualAux(2,i) = atoms%taual(2,i) - FLOOR(atoms%taual(2,i))
         taualAux(3,i) = atoms%taual(3,i) - FLOOR(atoms%taual(3,i))
         posAux(:,i) = MATMUL(amatAux,taualAux(:,i))
      END DO
      CALL inv3(amatAux,invAmatAux,amatAuxDet)



!     5. For the reference atom in auxiliary unit cell collect shortest distances
!        to other atoms in neighborhood

      ALLOCATE(sqrDistances(maxAtoms))
      ALLOCATE(neighborAtoms(maxAtoms))
      ALLOCATE(neighborAtomsDiff(3,maxAtoms))
      ALLOCATE(distIndexList(maxAtoms))
      ALLOCATE (nearestNeighbors(maxAtoms))
      ALLOCATE (nearestNeighborDists(maxAtoms))
      ALLOCATE (nearestNeighborDiffs(3,maxAtoms))

      !Find the reference atom
      iAtom = 0
      DO n = 1, atoms%ntype
         DO na = 1, atoms%neq(n)
            iAtom = iAtom + 1
            IF((n.EQ.refAtom).AND.na.EQ.1) THEN
               refPos(:) = posAux(:,iAtom)
            ENDIF
         ENDDO
      ENDDO

      !Collect the Distances between the refAtom and all other atoms
      neighborAtoms = 0
      iNeighborAtom = 0
      identicalAtoms = 0
      DO i = cubeStartIndex, cubeEndIndex
         DO j = cubeStartIndex, cubeEndIndex
            DO k = cubeStartIndex, cubeEndIndex
               DO m = 1, 3
                  offsetPos(m) = i*amatAux(m,1) + j*amatAux(m,2) + k*amatAux(m,3)
               END DO
               iAtom = 0
               DO n = 1, atoms%ntype
                  DO na = 1, atoms%neq(n)
                     iAtom = iAtom + 1
                     pos(:) = posAux(:,iAtom) + offsetPos(:)
                     currentDist = (refPos(1) - pos(1))**2 + &
                                   (refPos(2) - pos(2))**2 + &
                                   (refPos(3) - pos(3))**2
                     currentDiff = refPos(:) - pos(:)
                     IF (currentDist.LT.0.000001) THEN
                        identicalAtoms = identicalAtoms + 1
                     ELSE
                        iNeighborAtom = iNeighborAtom + 1
                        neighborAtoms(iNeighborAtom) = n
                        neighborAtomsDiff(:,iNeighborAtom) = currentDiff(:)
                        sqrDistances(iNeighborAtom) = currentDist
                     END IF
                  ENDDO
               END DO
            END DO
         END DO
      END DO
      IF (identicalAtoms.GT.1) THEN
         WRITE(*,*) 'Position: ', refPos(:)
         CALL juDFT_error("Too many atoms at same position.",calledby ="addNearestNeighbours_gfelem")
      END IF

      !Sort the atoms according to distance
      numNearestNeighbors = iNeighborAtom
      CALL sort(distIndexList(:iNeighborAtom),sqrDistances(:iNeighborAtom))
      DO i = 1, numNearestNeighbors
         nearestNeighbors(i) = neighborAtoms(distIndexList(i))
         nearestNeighborDists(i) = SQRT(sqrDistances(distIndexList(i)))
         nearestNeighborDiffs(:,i) = neighborAtomsDiff(:,distIndexList(i))
      END DO
      DEALLOCATE(sqrDistances,distIndexList,neighborAtomsDiff,neighborAtoms)

      !Maximum number of shells is number of atoms
      ALLOCATE(shellDistance(maxAtoms),source = 0.0)
      ALLOCATE(shellDiff(3,maxAtoms,maxAtoms),source = 0.0)
      ALLOCATE(shellAtom(maxAtoms),source=0)
      ALLOCATE(numshellAtoms(maxAtoms),source=0)

      !Sort the nearestNeighbours into shells
      lastIndex = 1 !Skip the first element (onsite)
      ishell = 1
      DO

         minDist = MINVAL(nearestNeighborDists(lastIndex:numNearestNeighbors))
         shellDistance(ishell) = minDist
         numshellAtoms(ishell) = 0
         DO iAtom = lastIndex, numNearestNeighbors
            lastIndex = iAtom
            IF(ABS(nearestNeighborDists(iAtom)-minDist).GT.1e-12) EXIT !List is sorted
            numshellAtoms(ishell) = numshellAtoms(ishell) + 1
            IF(shellAtom(ishell) == 0) THEN
               shellAtom(ishell) = nearestNeighbors(iAtom)
            ELSE IF(shellAtom(ishell) .NE. nearestNeighbors(iAtom)) THEN
               CALL juDFT_error("Found inequivalent atoms at same distance (not yet implemented)"&
                                ,calledby ="addNearestNeighbours_gfelem")
            ENDIF
            shellDiff(:,numshellAtoms(ishell),ishell) = nearestNeighborDiffs(:,iAtom)
         ENDDO

         IF (lastIndex<numNearestNeighbors) THEN
            ishell = ishell + 1
         ELSE
            EXIT
         ENDIF
      ENDDO
      DEALLOCATE(nearestNeighborDiffs,nearestNeighborDists,nearestNeighbors)


      ALLOCATE(shellAux(3,maxAtoms),source=0.0)
      ALLOCATE(shellAux1(3,maxAtoms),source=0.0)
      nshellsFound = ishell !We only want to consider nshells
      !Symmetry reduction (modernized and modified version of nshell.f from v26)
      DO ishell = 1, SIZE(shellDiff,3)
         IF(ishell.GT.nshellsFound) EXIT !We have finished the requested shells
         nshellDist = 0
         lastDist = 0.0
         DO ishell1 = 1, ishell
            IF(shellDistance(ishell1)-lastDist.GT.1e-12) nshellDist = nshellDist + 1
            lastDist = shellDistance(ishell1)
         ENDDO
         IF(nshellDist>nshells) EXIT
         !Take the representative element of the shell
         shellAux = 0.0
         shellAux(:,1) = shellDiff(:,1,ishell)

         nshellAtom = 1
         symLoop: DO iop = 1, sym%nop
            diffRot = matmul(sym%mrot(:,:,iop),shellAux(:,1))

            DO ishellAtom = 1, nshellAtom
               !Is the atom equivalent to another atom already in the shell
               IF(ALL(ABS(diffRot-shellAux(:,ishellAtom)).LT.tol)) CYCLE symLoop
            ENDDO

            nshellAtom = nshellAtom + 1
            shellAux(:,nshellAtom) = diffRot
         ENDDO symLoop

         IF(nshellAtom < numshellAtoms(ishell)) THEN  !Not all elements can be constructed from the representative element

            shellAux1 = 0.0
            !Find the atoms which are not represented
            nshellAtom1 = 0
            atomLoop: DO iAtom = 1, numshellAtoms(ishell)
               DO ishellAtom = 1, nshellAtom
                  IF(ALL(ABS(shellDiff(:,iAtom,ishell)-shellAux(:,ishellAtom)).LT.tol)) CYCLE atomLoop!.OR.&
                     !ALL(ABS(shellDiff(:,iAtom,ishell)+shellAux(:,ishellAtom)).LT.tol)) CYCLE atomLoop
               ENDDO

               nshellAtom1 = nshellAtom1 + 1
               shellAux1(:,nshellAtom1) = shellDiff(:,iAtom,ishell)
            ENDDO atomLoop

            !We have found a new shell
            nshellsFound = nshellsFound + 1
            IF(nshellsFound > SIZE(shellDiff,3)) CALL juDFT_error("Dim Error: nshells > SIZE",calledby ="addNearestNeighbours_gfelem")

            !Reorder shellDiff and all other arrays array (move everything above ishell one up to make space)
            DO ishell1 = nshellsFound, ishell + 2, -1
               shellAtom(ishell1) = shellAtom(ishell1-1)
               shellDistance(ishell1) = shellDistance(ishell1-1)
               numshellAtoms(ishell1) = numshellAtoms(ishell1-1)
               shellDiff(:,:,ishell1) = shellDiff(:,:,ishell1-1)
            ENDDO

            !Modify ishell (fewer atoms)
            numshellAtoms(ishell) = nshellAtom
            DO ishellAtom = 1, nshellAtom
               shellDiff(:,ishellAtom,ishell) = shellAux(:,ishellAtom)
            ENDDO

            !Insert Element at ishell+1 (This way it will be the next element in the
            !loop if it needs to be deconstructed further)
            shellAtom(ishell+1) = shellAtom(ishell)
            shellDistance(ishell+1) = shellDistance(ishell)
            numshellAtoms(ishell+1) = nshellAtom1
            DO ishellAtom = 1, nshellAtom1
               shellDiff(:,ishellAtom,ishell+1) = shellAux1(:,ishellAtom)
            ENDDO

         ENDIF

      ENDDO

      nshellsFound = ishell - 1

      DO ishell = 1, nshellsFound
         IF(l_write) THEN
            WRITE(oUnit,'(/,A,f14.8)') 'Adding shell with distance: ', shellDistance(ishell)

            WRITE(oUnit,'(/,A)') ' Contains the following atom pairs:'
            DO ishellAtom = 1, numshellAtoms(ishell)
               WRITE(oUnit,'(3f14.8)') shellDiff(:,ishellAtom,ishell)
            ENDDO
         ENDIF

         !Transform represenative element to lattice coordinates
         diff = MATMUL(invAmatAux,shellDiff(:,1,ishell))
         !l_sphavg has to be false
         i_gf =  this%add(l,refAtom,iContour,.FALSE.,lp=lp,nTypep=shellAtom(ishell),&
                          atomDiff=diff,l_fixedCutoffset=l_fixedCutoffset,&
                          fixedCutoff=fixedCutoff)

         this%elem(i_gf)%refCutoff = refCutoff

         IF(l_write) THEN
            WRITE(oUnit,'(A,I6,I6,6f14.8)') 'GF Element: ', refAtom, shellAtom(ishell),&
                                            shellDiff(:,1,ishell), diff(:)
         ENDIF

         !Add negative of diff (for Jij we need Gij and Gji (could maybe be done with conjugation))
         !This should not produce problems if symmertry is reduced because add makes sure that there
         !are no duplicates in this%elem
         diff = -1.0 * diff
         i_gf =  this%add(l,shellAtom(ishell),iContour,.FALSE.,lp=lp,nTypep=refAtom,&
                          atomDiff=diff,l_fixedCutoffset=l_fixedCutoffset,&
                          fixedCutoff=fixedCutoff)

         this%elem(i_gf)%refCutoff = refCutoff

         IF(l_write) THEN
            WRITE(oUnit,'(A,I6,I6,6f14.8)') 'GF Element: ', shellAtom(ishell), refAtom, &
                                            -1.0 * shellDiff(:,1,ishell), diff(:)
         ENDIF

      ENDDO

   CALL timestop("Green's Function: Add nearest Neighbors")

   END SUBROUTINE addNearestNeighbours_gfelem

   INTEGER FUNCTION find_gfelem(this,l,nType,iContour,l_sphavg,lp,nTypep,atomDiff,uniqueMax,l_found) result(i_gf)

      !Maps between the four indices (l,lp,nType,nTypep) and the position in the
      !gf arrays

      CLASS(t_gfinp),      INTENT(IN)    :: this
      INTEGER,             INTENT(IN)    :: l
      INTEGER,             INTENT(IN)    :: nType
      INTEGER,             INTENT(IN)    :: iContour
      LOGICAL,             INTENT(IN)    :: l_sphavg
      INTEGER, OPTIONAL,   INTENT(IN)    :: lp
      INTEGER, OPTIONAL,   INTENT(IN)    :: nTypep
      REAL,    OPTIONAL,   INTENT(IN)    :: atomDiff(:)
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
            IF(nTypep/=0) THEN
               IF(this%elem(i_gf)%atomTypep.NE.nTypep) CYCLE
            ELSE
               IF(this%elem(i_gf)%atomTypep.NE.nType) CYCLE
            ENDIF
         ELSE
            IF(this%elem(i_gf)%atomTypep.NE.nType) CYCLE
         ENDIF
         IF(this%elem(i_gf)%l_sphavg .neqv. l_sphavg) CYCLE
         !Check the phasefactor
         IF(PRESENT(atomDiff)) THEN
            IF(ABS(this%elem(i_gf)%atomDiff(1)-atomDiff(1)).GT.1e-12.OR.&
               ABS(this%elem(i_gf)%atomDiff(2)-atomDiff(2)).GT.1e-12.OR.&
               ABS(this%elem(i_gf)%atomDiff(3)-atomDiff(3)).GT.1e-12) CYCLE
         ELSE
            IF(ABS(this%elem(i_gf)%atomDiff(1)).GT.1e-12.OR.&
               ABS(this%elem(i_gf)%atomDiff(2)).GT.1e-12.OR.&
               ABS(this%elem(i_gf)%atomDiff(3)).GT.1e-12) CYCLE
         ENDIF
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

   PURE LOGICAL FUNCTION checkOnsite_gfinp(this)

      !Check if there are any oniste elements
      CLASS(t_gfinp),               INTENT(IN)    :: this

      INTEGER :: i_gf

      checkOnsite_gfinp = .FALSE.
      DO i_gf = 1, this%n
         IF(this%elem(i_gf)%l.NE.this%elem(i_gf)%lp) CYCLE
         IF(this%elem(i_gf)%atomType.NE.this%elem(i_gf)%atomTypep) CYCLE
         IF(ANY(ABS(this%elem(i_gf)%atomDiff).GT.1e-12)) CYCLE
         checkOnsite_gfinp = .TRUE.
      ENDDO
   END FUNCTION checkOnsite_gfinp

   PURE LOGICAL FUNCTION checkOffdiagonal_gfinp(this)

      !Check if there are any oniste elements
      CLASS(t_gfinp),               INTENT(IN)    :: this

      INTEGER :: i_gf

      checkOffdiagonal_gfinp = .FALSE.
      DO i_gf = 1, this%n
         IF(this%elem(i_gf)%l.EQ.this%elem(i_gf)%lp.AND. &
            this%elem(i_gf)%atomType.EQ.this%elem(i_gf)%atomTypep.AND. &
            ALL(ABS(this%elem(i_gf)%atomDiff).LT.1e-12)) CYCLE
         checkOffdiagonal_gfinp = .TRUE.
      ENDDO
   END FUNCTION checkOffdiagonal_gfinp

   PURE INTEGER FUNCTION countLOs_gfelem(this,atoms)

      !Counts the number of LOs associated with this green's function element
      USE m_types_atoms

      CLASS(t_gfelementtype),   INTENT(IN)  :: this
      TYPE(t_atoms),            INTENT(IN)  :: atoms

      INTEGER :: l,lp,atomType,atomTypep,ilo,ilop

      l  = this%l
      lp = this%lp
      atomType  = this%atomType
      atomTypep = this%atomTypep

      countLOs_gfelem = 0
      DO ilo = 1, atoms%nlo(atomType)
         IF(atoms%llo(ilo,atomType).NE.l) CYCLE
         countLOs_gfelem = countLOs_gfelem + 1
      ENDDO

      IF(l.NE.lp.OR.atomType.NE.atomTypep) THEN
         DO ilop = 1, atoms%nlo(atomTypep)
            IF(atoms%llo(ilop,atomType).NE.lp) CYCLE
            countLOs_gfelem = countLOs_gfelem + 1
         ENDDO
      ENDIF

   END FUNCTION countLOs_gfelem

   PURE LOGICAL FUNCTION isOffDiag_gfelem(this)

      CLASS(t_gfelementtype),   INTENT(IN)  :: this

      isoffDiag_gfelem = this%l.NE.this%lp.OR.this%atomType.NE.this%atomTypep&
                           .OR.ANY(ABS(this%atomDiff).GT.1e-12)

   END FUNCTION isOffDiag_gfelem

END MODULE m_types_gfinp
