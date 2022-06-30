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

   REAL,    PARAMETER :: ATOMDIFF_EPS = 1e-5

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

      !Symmetry relations to other gf element (only intersite)
      INTEGER :: representative_elem = -1
      INTEGER :: representative_op = -1
      REAL    :: representative_diff(3)  = [0.0,0.0,0.0] !Distance between atoms of representative element (lattice coordinates) for intersite phase
      INTEGER :: atom = 0   !Specific atom for this element (Used for intersite elements)
      INTEGER :: atomp = 0  !Specific atom for this element (Used for intersite elements)

      !K-resolved switches
      LOGICAL :: l_kresolved = .FALSE. !Should the Greens function be calculated k-resolved
      LOGICAL :: l_kresolved_int = .FALSE. !Should the Greens function be calculated k-resolved up after the Kramers-Kronig
                                           !Transformation (Intersite elements)
   CONTAINS
      PROCEDURE :: init                => init_gfelem
      PROCEDURE :: countLOs            => countLOs_gfelem !Count the local orbitals attached to the element
      PROCEDURE :: isOffDiag           => isOffDiag_gfelem !Is this element offdiagonal (i.e either l/=lp or intersite)
      PROCEDURE :: isIntersite         => isIntersite_gfelem !Is this element intersite (either unequal atomypes or non zero atomDiff)
      PROCEDURE :: equals              => equals_gfelem !Is the element equal to another (For deduplicating added elements)
      PROCEDURE :: equals_coefficients => equals_coefficients_gfelem !Is the element equal to another from the perspective of the BZ Coefficients
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
      INTEGER, ALLOCATABLE :: torqueElem(:,:)
      INTEGER, ALLOCATABLE :: numTorqueElems(:)
      INTEGER, ALLOCATABLE :: intersiteAtomicNumberSelection(:,:)
      INTEGER, ALLOCATABLE :: numintersiteAtomicNumberSelection(:)
   CONTAINS
      PROCEDURE :: read_xml             => read_xml_gfinp
      PROCEDURE :: mpi_bc               => mpi_bc_gfinp
      PROCEDURE :: distribute_elements  => distribute_elements_gfinp
      PROCEDURE :: init                 => init_gfinp
      PROCEDURE :: find_gfelem_simple
      PROCEDURE :: find_gfelem_type
      GENERIC   :: find                 => find_gfelem_simple, find_gfelem_type
      PROCEDURE :: find_symmetry_rotated_bzcoeffs => find_symmetry_rotated_bzcoeffs_gfinp
      PROCEDURE :: find_contour         => find_contour
      PROCEDURE :: add                  => add_gfelem
      PROCEDURE :: addNearestNeighbours => addNearestNeighbours_gfelem
      PROCEDURE :: getuniqueElement     => getuniqueElement_gfinp
      PROCEDURE :: isUnique             => isUnique_gfinp
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
      CALL mpi_bc(this%minCalcDistance,rank,mpi_comm)
      CALL mpi_bc(this%l_outputSphavg,rank,mpi_comm)
      CALL mpi_bc(this%l_intFullRadial,rank,mpi_comm)
      CALL mpi_bc(this%n,rank,mpi_comm)
      CALL mpi_bc(this%ne,rank,mpi_comm)
      CALL mpi_bc(this%ellow,rank,mpi_comm)
      CALL mpi_bc(this%elup,rank,mpi_comm)
      CALL mpi_bc(this%numberContours,rank,mpi_comm)
      CALL mpi_bc(this%hiaElem,rank,mpi_comm)
      CALL mpi_bc(this%torqueElem,rank,mpi_comm)
      CALL mpi_bc(this%numTorqueElems,rank,mpi_comm)
      call mpi_bc(this%intersiteAtomicNumberSelection,rank,mpi_comm)
      call mpi_bc(this%numintersiteAtomicNumberSelection,rank,mpi_comm)

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
         CALL mpi_bc(this%elem(n)%representative_elem,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%representative_op,rank,mpi_comm)
         CALL mpi_bc(rank,mpi_comm,this%elem(n)%representative_diff)
         CALL mpi_bc(this%elem(n)%atom,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%atomp,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%l_kresolved,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%l_kresolved_int,rank,mpi_comm)
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

   SUBROUTINE distribute_elements_gfinp(this, rank, size, nspins, i_gf_start, i_gf_end, spin_start, spin_end, k_resolved)
      !Distribute the Greens function elements that are not kresolved for the Kramers Kronig
      !integration
      CLASS(t_gfinp),    INTENT(IN)  :: this
      INTEGER,           INTENT(IN)  :: rank, size, nspins
      INTEGER,           INTENT(OUT) :: i_gf_start, i_gf_end, spin_start, spin_end
      LOGICAL, OPTIONAL, INTENT(IN)  :: k_resolved

      INTEGER :: n_elems, i_gf, currentIndex, n_gf_task, extra
      INTEGER :: i_elem_start, i_elem_end
      LOGICAL :: k_resolved_arg

      k_resolved_arg = .FALSE.
      IF(PRESENT(k_resolved)) k_resolved_arg = k_resolved


      n_elems = COUNT(.NOT.this%elem(:)%l_kresolved_int)
      IF(k_resolved_arg) THEN
         n_elems = this%n - n_elems
      ENDIF

#ifdef CPP_MPI
      IF(size>1) THEN
         IF(n_elems>=size) THEN
            !Just distribute the individual gf elements over the ranks
            n_gf_task = FLOOR(REAL(n_elems)/(size))
            extra = n_elems - n_gf_task*size
            i_elem_start = rank*n_gf_task + 1 + extra
            i_elem_end = (rank+1)*n_gf_task   + extra
            IF(rank < extra) THEN
               i_elem_start = i_elem_start - (extra - rank)
               i_elem_end = i_elem_end - (extra - rank - 1)
            ENDIF
            spin_start = 1
            spin_end   = nspins
         ELSE IF(n_elems*nspins>size) THEN
            !Just fill up the ranks
            i_elem_start = rank + 1
            i_elem_end   = rank + 1
            spin_start = 1
            spin_end   = nspins
         ELSE
            !If there are few enough gf elements then distribute the spins
            spin_start = MOD(rank,nspins) + 1
            spin_end   = MOD(rank,nspins) + 1
            i_elem_start = 1 + FLOOR(REAL(rank)/nspins)
            i_elem_end   = 1 + FLOOR(REAL(rank)/nspins)
         ENDIF
      ELSE
         !Distribute nothing
         i_elem_start = 1
         i_elem_end = n_elems
         spin_start = 1
         spin_end   = nspins
      ENDIF
#else
      i_elem_start = 1
      i_elem_end = n_elems
      spin_start = 1
      spin_end   = nspins
#endif

      IF(i_elem_start.LT.1 .OR. i_elem_start.GT.n_elems.OR.&
         i_elem_end.LT.1 .OR. i_elem_end.GT.n_elems) THEN
         i_gf_start = -1
         i_gf_end = -1
         RETURN
      ENDIF

      currentIndex = 0
      DO i_gf = 1, this%n
         IF(k_resolved_arg) THEN
            IF(.NOT.this%elem(i_gf)%l_kresolved_int) CYCLE
         ELSE
            IF(this%elem(i_gf)%l_kresolved_int) CYCLE
         ENDIF
         currentIndex = currentIndex + 1
         IF(currentIndex.EQ.i_elem_start) THEN
            i_gf_start = i_gf
         ENDIF
         IF(currentIndex.EQ.i_elem_end) THEN
            i_gf_end = i_gf
            EXIT
         ELSE IF(i_gf.EQ.this%n) THEN
            CALL juDFT_error('Distribution of Greens functions elements failed', calledby='distribute_elements_gfinp')
         ENDIF
      ENDDO

   END SUBROUTINE distribute_elements_gfinp

   SUBROUTINE read_xml_gfinp(this, xml)
      USE m_types_xml
      CLASS(t_gfinp), INTENT(INOUT):: this
      TYPE(t_xml),INTENT(INOUT) ::xml

      INTEGER :: numberNodes,ntype,itype,n_hia,i_gf,refL,refGF,nshells
      INTEGER :: i,l,lp,iContour,iContourp, atomicNumber, shellAtomicNumberSelection
      REAL    :: fixedCutoff
      CHARACTER(len=200)  :: xPathA,xPathS,label,cutoffArg,str, shellElement
      CHARACTER(len=1),PARAMETER :: spdf(0:3) = ['s','p','d','f']
      LOGICAL :: l_gfinfo_given,l_fixedCutoffset,l_sphavg,l_kresolved,l_kresolved_int, l_found
      LOGICAL :: lp_calc(0:3,0:3)

      xPathA = '/fleurInput/calculationSetup/greensFunction'
      numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
      l_gfinfo_given = numberNodes.EQ.1

      IF (l_gfinfo_given) THEN
         this%l_mperp=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
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

      ALLOCATE(this%hiaElem(4*ntype))
      ALLOCATE(this%numTorqueElems(ntype),source=0)
      ALLOCATE(this%torqueElem(ntype,(lmaxU_const+1)**2),source=-1)
      ALLOCATE(this%intersiteAtomicNumberSelection(ntype,(lmaxU_const+1)**2),source=-1)
      ALLOCATE(this%numintersiteAtomicNumberSelection((lmaxU_const+1)**2*ntype),source=0)

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
            IF(xml%versionNumber>=35) THEN
               l_kresolved = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@k_resolved'))
            ELSE
               l_kresolved = .FALSE.
            ENDIF
            IF(xml%versionNumber>=36) THEN
               shellElement = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@shellElement')))
            ELSE
               shellElement = 'all'
            ENDIF
            shellAtomicNumberSelection = -1
            if (trim(adjustl(shellElement))/='all') then
               do atomicNumber = 1, size(namat_const)
                  if (trim(adjustl(shellElement))==namat_const(atomicNumber)) then
                     shellAtomicNumberSelection = atomicNumber
                     exit
                  else if(atomicNumber == size(namat_const)) then
                     CALL juDFT_error("Error reading in gf-information: Wrong Atomic number selection given value is not a chemical element",&
                                      calledby="read_xml_gfinp")
                  endif
               enddo
            endif
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
                                      fixedCutoff=fixedCutoff,nshells=nshells,k_resolved=l_kresolved)
                     if (nshells /= 0 .and. shellAtomicNumberSelection /= -1) THEN
                        this%numintersiteAtomicNumberSelection(i_gf) = this%numintersiteAtomicNumberSelection(i_gf) + 1
                        this%intersiteAtomicNumberSelection(i_gf,this%numintersiteAtomicNumberSelection(i_gf)) = shellAtomicNumberSelection
                     endif                 
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
                                   fixedCutoff=fixedCutoff,nshells=nshells,k_resolved=l_kresolved)
                  if (nshells /= 0 .and. shellAtomicNumberSelection /= -1) THEN
                     this%numintersiteAtomicNumberSelection(i_gf) = this%numintersiteAtomicNumberSelection(i_gf) + 1
                     this%intersiteAtomicNumberSelection(i_gf,this%numintersiteAtomicNumberSelection(i_gf)) = shellAtomicNumberSelection
                  endif
               ENDDO
            ENDIF

            !Find the reference element
            IF(refL /= -1 .AND.ANY(lp_calc)) THEN
               !Find the element
               refGF = this%find(refL,itype,iContour,l_sphavg,k_resolved=.FALSE.,k_resolved_int=.FALSE.,nshells=nshells,l_found=l_found)
               IF(.NOT.l_found) THEN
                  refGF = this%add(refL,itype,iContour,l_sphavg,k_resolved=.FALSE.)
               ENDIF
               DO l = 0,lmaxU_const
                  DO lp = 0,lmaxU_const
                     IF(.NOT.lp_calc(lp,l)) CYCLE
                     i_gf = this%find(l,itype,iContour,l_sphavg,lp=lp,nshells=nshells,k_resolved=l_kresolved)
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
                             fixedCutoff=fixedCutoff,k_resolved=.FALSE.)
            n_hia = n_hia + 1
            this%hiaElem(n_hia) = i_gf
         ENDDO

         WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/torqueCalculation'
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

            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/torqueCalculation/greensfElements'
            DO l = 0,lmaxU_const
               str = xml%GetAttributeValue(TRIM(xPathA)//'/'//spdf(l))
               READ(str,'(4l2)') (lp_calc(lp,l),lp=0,3)
               DO lp = 0,lmaxU_const
                  IF(.NOT.lp_calc(lp,l)) CYCLE
                  !Torque GF has to have radial dependence
                  i_gf =  this%add(l,itype,iContour,.FALSE.,lp=lp,l_fixedCutoffset=l_fixedCutoffset,&
                                   fixedCutoff=fixedCutoff,k_resolved=.FALSE.)
                  this%numTorqueElems(itype) = this%numTorqueElems(itype) + 1
                  this%torqueElem(itype,this%numTorqueElems(itype)) = i_gf
               ENDDO
            ENDDO

            !Find the reference element
            IF(refL /= -1 .AND.ANY(lp_calc)) THEN
               !Find the element
               refGF = this%find(refL,itype,iContour,.FALSE.,k_resolved=.FALSE.,k_resolved_int=.FALSE., l_found=l_found)
               IF(.NOT.l_found) THEN
                  refGF = this%add(refL,itype,iContour,.FALSE.,k_resolved=.FALSE.)
               ENDIF
               DO l = 0,lmaxU_const
                  DO lp = 0,lmaxU_const
                     IF(.NOT.lp_calc(lp,l)) CYCLE
                     i_gf = this%find(l,itype,iContour,.FALSE.,lp=lp,k_resolved=.FALSE.,k_resolved_int=.FALSE.)
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

      INTEGER :: i_gf,l,lp,atomType,iContour,refCutoff
      INTEGER :: refCutoff1,nOtherAtoms,nOtherAtoms1,iOtherAtom,lref,n_intersite, i_inter
      LOGICAL :: l_sphavg, l_all_kresolved, l_kresolved_radial
      INTEGER :: hiaElem(atoms%n_hia), intersite_elems(this%n), shells(this%n)
      LOGICAL :: written(atoms%nType), l_kresolved
      TYPE(t_gfelementtype), ALLOCATABLE :: gfelem(:)
      INTEGER, ALLOCATABLE :: atomTypepList(:),atomTypepList1(:), atomicNumberSelection(:)

      IF(this%n==0) RETURN !Nothing to do here

      n_intersite = 0
      shells = 0
      intersite_elems = -1
      written = .FALSE.
      !Find the elements for which we need to compute the nearest neighbours
      DO i_gf = 1, this%n
         IF(this%elem(i_gf)%atomTypep<0) THEN
            n_intersite = n_intersite + 1
            intersite_elems(n_intersite) = i_gf
            shells(n_intersite) = ABS(this%elem(i_gf)%atomTypep)
            this%elem(i_gf)%atomTypep = this%elem(i_gf)%atomType
         ENDIF
      ENDDO
      DO i_inter = 1, n_intersite
         i_gf = intersite_elems(i_inter)
         l  = this%elem(i_gf)%l
         lp = this%elem(i_gf)%lp
         atomType  = this%elem(i_gf)%atomType
         iContour = this%elem(i_gf)%iContour
         refCutoff = this%elem(i_gf)%refCutoff
         l_sphavg = this%elem(i_gf)%l_sphavg
         l_kresolved = this%elem(i_gf)%l_kresolved

         refCutoff = MERGE(i_gf,refCutoff,refCutoff==-1) !If no refCutoff is set for the intersite element
                                                         !we take the onsite element as reference
         atomicNumberSelection = this%intersiteAtomicNumberSelection(i_gf,:this%numintersiteAtomicNumberSelection(i_gf))
         CALL this%addNearestNeighbours(shells(i_inter),l,lp,atomType,l_sphavg,iContour,l_kresolved,this%elem(i_gf)%l_fixedCutoffset,&
                                        this%elem(i_gf)%fixedCutoff,refCutoff,atoms,cell,sym,input,&
                                        .NOT.written(atomType),nOtherAtoms,atomTypepList,&
                                        atomicNumberSelection=atomicNumberSelection)
         written(atomType) = .TRUE.

         !Add the other atomtypes (i,j) -> (j,i)
         DO iOtherAtom = 1, nOtherAtoms
            atomType = atomTypepList(iOtherAtom)
            !First add the reference cutoff element
            lref = this%elem(refCutoff)%l
            refCutoff1 =  this%add(lref,atomType,iContour,l_sphavg,k_resolved=.FALSE.,&
                                   l_fixedCutoffset=this%elem(i_gf)%l_fixedCutoffset,&
                                   fixedCutoff=this%elem(i_gf)%fixedCutoff)

            WRITE(oUnit,'(A,i0)') 'Adding shells for atom: ', atomType
            
            if (size(atomicNumberSelection) /= 0) then
               CALL this%addNearestNeighbours(shells(i_inter),l,lp,atomType,l_sphavg,iContour,l_kresolved,this%elem(i_gf)%l_fixedCutoffset,&
                                             this%elem(i_gf)%fixedCutoff,refCutoff1,atoms,cell,sym,input,&
                                             .NOT.written(atomType),nOtherAtoms1,atomTypepList1,&
                                             atomicNumberSelection=(/atoms%nz(atomType)/))
            else
               CALL this%addNearestNeighbours(shells(i_inter),l,lp,atomType,l_sphavg,iContour,l_kresolved,this%elem(i_gf)%l_fixedCutoffset,&
                                             this%elem(i_gf)%fixedCutoff,refCutoff1,atoms,cell,sym,input,&
                                             .NOT.written(atomType),nOtherAtoms1,atomTypepList1)
            endif

         ENDDO
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

      l_all_kresolved = .FALSE.
      l_kresolved_radial = .FALSE.
      DO i_gf = 1, this%n
         l_sphavg  = this%elem(i_gf)%l_sphavg
         IF(this%elem(i_gf)%l_kresolved) THEN
            l_all_kresolved = .TRUE.
            IF(.NOT.l_sphavg) THEN
               l_kresolved_radial = .TRUE.
            ENDIF
         ENDIF

      ENDDO

      IF(this%minCalcDistance>=0.0) THEN
         IF(input%mindistance>this%minCalcDistance) THEN
            CALL juDFT_warn("The minimum Distance for Green's Function Calculation"// &
                            "is smaller than the distance requirement:"//&
                            "No Green's Functions will be calculated", calledby="init_gfinp")
         ENDIF
      ENDIF

      IF(l_all_kresolved.AND.input%bz_integration/=BZINT_METHOD_HIST) THEN
         CALL juDFT_error("Completely k-resolved Greens functions only implemented for histogram method",&
                          calledby="init_gfinp")
      ENDIF

      IF(l_kresolved_radial) THEN
         CALL juDFT_error("k-resolved Greens functions + Radial dependence not implemented",&
                          calledby="init_gfinp")
      ENDIF

      IF(ANY(this%numTorqueElems(:)>0)) THEN
         IF(input%jspins.NE.2) CALL juDFT_error("Torque calculation only for magnetic systems", calledby="init_gfinp")
         IF(.NOT.noco%l_mperp.or.this%l_mperp) &
            CALL juDFT_error("Torque calculation only with l_mperp=T (both noco and greensfunction)", &
                             calledby="init_gfinp")
         IF(sym%nop>1) CALL juDFT_warn("Torque calculation only without symmetries", calledby="init_gfinp")
      ENDIF


      WRITE(oUnit,'(/,A,I0)') "Green's Function Elements: ", this%n
      WRITE(oUnit,'(A)') "Index | l/lp | atom/atomp | contour | sphavg | refCutoff | repr_elem(repr_op) | k_resolved | atomDiff"
      WRITE(oUnit,'(A)') "-----------------------------------------------------------------------------------------------------------------"
      DO i_gf = 1, this%n
         WRITE(oUnit,9000) i_gf, this%elem(i_gf)%l,this%elem(i_gf)%lp,this%elem(i_gf)%atomType,this%elem(i_gf)%atomTypep,&
                           this%elem(i_gf)%iContour,this%elem(i_gf)%l_sphavg,this%elem(i_gf)%refCutoff,&
                           this%elem(i_gf)%representative_elem,this%elem(i_gf)%representative_op, &
                           this%elem(i_gf)%l_kresolved, this%elem(i_gf)%l_kresolved_int, &
                           this%elem(i_gf)%atomDiff(:)
      ENDDO
      WRITE(oUnit,'(/)')
9000  FORMAT(I5, " | ",I1,"/",I1,"  |",I5,"/",I5," | ",I7," | ",l6," | ", I9," | ", I9,"(",I2,")      | ",l6,"(",l1,")  |",3f7.3)

   END SUBROUTINE init_gfinp

   PURE LOGICAL FUNCTION isUnique_gfinp(this,index, distinct_kresolved_int, distinct_symmetry_equivalent_diffs)
      !Return whether the given element is the first with the combination
      !of l lp, atomType, atomTypep, l_sphavg, l_kresolved

      CLASS(t_gfinp),   INTENT(IN)  :: this
      INTEGER,          INTENT(IN)  :: index
      LOGICAL, OPTIONAL,INTENT(IN)  :: distinct_kresolved_int
      LOGICAL, OPTIONAL,INTENT(IN)  :: distinct_symmetry_equivalent_diffs

      INTEGER :: i_gf, uniqueIndex

      uniqueIndex = this%getuniqueElement(index, distinct_kresolved_int, distinct_symmetry_equivalent_diffs)
      isunique_gfinp = uniqueIndex == index

   END FUNCTION isUnique_gfinp

   PURE INTEGER FUNCTION getuniqueElement_gfinp(this, index, distinct_kresolved_int, distinct_symmetry_equivalent_diffs) Result(uniqueIndex)

      CLASS(t_gfinp),   INTENT(IN)  :: this
      INTEGER,          INTENT(IN)  :: index
      LOGICAL, OPTIONAL,INTENT(IN)  :: distinct_kresolved_int
      LOGICAL, OPTIONAL,INTENT(IN)  :: distinct_symmetry_equivalent_diffs

      DO uniqueIndex = 1, index
         !If the element has a representative element set it can not be unique
         !IF(this%elem(uniqueIndex)%representative_elem>0) CYCLE
         IF(this%elem(uniqueIndex)%equals_coefficients(this%elem(index), distinct_kresolved_int, distinct_symmetry_equivalent_diffs)) THEN
            RETURN
         ENDIF
      ENDDO

   END FUNCTION getuniqueElement_gfinp

   INTEGER FUNCTION uniqueElements_gfinp(this,atoms, max_index, l_sphavg, lo, l_kresolved_int,maxLO) Result(uniqueElements)

      USE m_types_atoms

      CLASS(t_gfinp),   INTENT(IN)     :: this
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      INTEGER, OPTIONAL,INTENT(IN)     :: max_index
      LOGICAL, OPTIONAL,INTENT(IN)     :: l_sphavg !uniqueElements are determined separately for radial dependence and spherically averaging
      LOGICAL, OPTIONAL,INTENT(IN)     :: lo       !We are interested in unique LO elems (radial dependence)
      LOGICAL, OPTIONAL,INTENT(IN)     :: l_kresolved_int !K-resolved until the Kramers Kronig Integration
      INTEGER, OPTIONAL,INTENT(INOUT)  :: maxLO    !Maximum number of Elements associated with a GF element

      INTEGER :: maxGF,nLO
      INTEGER :: i_gf
      LOGICAL :: l_sphavgArg, l_sphavgElem,loArg, resolvedArg, l_kresolved_int_elem
      LOGICAL :: distinct_kresolved_int

      distinct_kresolved_int = PRESENT(l_kresolved_int)

      !Process optional switches and arguments
      l_sphavgArg = .TRUE.
      IF(PRESENT(l_sphavg)) l_sphavgArg = l_sphavg

      loArg = .FALSE.
      IF(PRESENT(lo)) loArg = lo

      uniqueElements = 0
      IF(PRESENT(maxLO)) maxLO = 0

      IF(PRESENT(max_index)) THEN
         maxGF = max_index
      ELSE
         maxGF = this%n
      ENDIF

      !Count the unique Elements before maxGF
      DO i_gf = 1, maxGF
         l_sphavgElem  = this%elem(i_gf)%l_sphavg
         l_kresolved_int_elem  = this%elem(i_gf)%l_kresolved_int

         IF(l_sphavgElem .neqv. l_sphavgArg) CYCLE

         IF(distinct_kresolved_int) THEN
            IF(l_kresolved_int_elem .neqv. l_kresolved_int) CYCLE
         ENDIF

         IF(.NOT.this%isUnique(i_gf,distinct_kresolved_int)) CYCLE

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

      ENDDO

   END FUNCTION uniqueElements_gfinp

   INTEGER FUNCTION add_gfelem(this,l,atomType,iContour,l_sphavg,lp,atomTypep,atomDiff,&
                               l_fixedCutoffset,fixedCutoff,nshells,k_resolved,atom,atomp) Result(i_gf)

      CLASS(t_gfinp),      INTENT(INOUT)  :: this
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: atomType
      INTEGER,             INTENT(IN)     :: iContour
      LOGICAL,             INTENT(IN)     :: l_sphavg
      INTEGER, OPTIONAL,   INTENT(IN)     :: lp
      INTEGER, OPTIONAL,   INTENT(IN)     :: atomTypep !Specify the second atom
      REAL,    OPTIONAL,   INTENT(IN)     :: atomDiff(:)
      LOGICAL, OPTIONAL,   INTENT(IN)     :: l_fixedCutoffset
      REAL,    OPTIONAL,   INTENT(IN)     :: fixedCutoff
      INTEGER, OPTIONAL,   INTENT(IN)     :: nshells
      LOGICAL, OPTIONAL,   INTENT(IN)     :: k_resolved
      INTEGER, OPTIONAL,   INTENT(IN)     :: atom, atomp

      LOGICAL l_found
      TYPE(t_gfelementtype) :: new_element
      TYPE(t_gfelementtype), ALLOCATABLE :: gfelem(:)

      CALL new_element%init(l,atomType,iContour,l_sphavg,lp=lp,atomTypep=atomTypep,&
                            nshells=nshells,atomDiff=atomDiff,k_resolved=k_resolved,&
                            l_fixedCutoffset=l_fixedCutoffset,fixedCutoff=fixedCutoff,&
                            atom=atom,atomp=atomp)

      !Check if this job has already been added
      i_gf = this%find(new_element,l_found=l_found)
      IF(l_found) RETURN

      IF(.NOT.ALLOCATED(this%elem)) THEN
         ALLOCATE(this%elem((lmaxU_const+1)**2))
      ENDIF

      this%n = this%n + 1

      IF(this%n>SIZE(this%elem)) THEN
         !Reallocate with doubled size
         ALLOCATE(gfelem(2*this%n))
         gfelem(:this%n-1) = this%elem(:this%n-1)
         CALL move_alloc(gfelem,this%elem)
      ENDIF

      this%elem(this%n) = new_element
      i_gf = this%n

   END FUNCTION add_gfelem

   INTEGER FUNCTION find_symmetry_rotated_bzcoeffs_gfinp(this, atoms, sym, i_gf, iop, l_sphavg, lo) RESULT(i_elem_rot)

      USE m_types_sym
      USE m_types_atoms

      CLASS(t_gfinp),         INTENT(IN)  :: this
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_sym),            INTENT(IN)  :: sym
      INTEGER,                INTENT(IN)  :: i_gf, iop
      LOGICAL,                INTENT(IN)  :: l_sphavg
      LOGICAL, OPTIONAL,      INTENT(IN)  :: lo

      TYPE(t_gfelementtype) :: gfelem_rot
      LOGICAL :: loArg
      REAL    :: diff(3)
      INTEGER :: atom_rot, atom_rotp, i_gf_rot, iop_arg


      IF(.NOT.this%elem(i_gf)%isIntersite()) CALL juDFT_error("find_symmetry_rotated_bzcoeffs should only be used for Intersite Green's functions", calledby='find_symmetry_rotated_bzcoeffs')

      iop_arg = iop
      IF(iop_arg > sym%nop) iop_arg = iop_arg - sym%nop

      gfelem_rot = this%elem(i_gf)

      diff = matmul(sym%mrot(:,:,iop_arg),this%elem(i_gf)%atomDiff)
      atom_rot = sym%mapped_atom(iop_arg, this%elem(i_gf)%atom)
      atom_rotp = sym%mapped_atom(iop_arg, this%elem(i_gf)%atomp)

      gfelem_rot%atomDiff = diff
      gfelem_rot%atom = atom_rot
      gfelem_rot%atomp = atom_rotp

      i_gf_rot = this%find(gfelem_rot,distinct_kresolved_int=.FALSE.)
      i_elem_rot = this%uniqueElements(atoms,max_index=i_gf_rot,l_sphavg=l_sphavg, lo=lo)

   END FUNCTION find_symmetry_rotated_bzcoeffs_gfinp

   SUBROUTINE addNearestNeighbours_gfelem(this,nshells,l,lp,refAtom,l_sphavg,iContour,l_kresolved,l_fixedCutoffset,fixedCutoff,&
                                          refCutoff,atoms,cell,sym,input,l_write,nOtherAtoms,atomTypepList,atomicNumberSelection)

      USE m_types_atoms
      USE m_types_cell
      USE m_types_sym
      USE m_types_input
      USE m_atom_shells

      !This is essentially a simplified version of chkmt, because we have a given
      !reference atom and do not need to consider all distances between all atoms

      CLASS(t_gfinp),      INTENT(INOUT)  :: this
      INTEGER,             INTENT(IN)     :: nshells !How many nearest neighbour shells are requested
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: lp
      INTEGER,             INTENT(IN)     :: refAtom !which is the reference atom
      LOGICAL,             INTENT(IN)     :: l_sphavg
      INTEGER,             INTENT(IN)     :: iContour
      LOGICAL,             INTENT(IN)     :: l_kresolved
      LOGICAL,             INTENT(IN)     :: l_fixedCutoffset
      REAL,                INTENT(IN)     :: fixedCutoff
      INTEGER,             INTENT(IN)     :: refCutoff
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_cell),        INTENT(IN)     :: cell
      TYPE(t_sym),         INTENT(IN)     :: sym
      TYPE(t_input),       INTENT(IN)     :: input
      LOGICAL,             INTENT(IN)     :: l_write
      INTEGER,             INTENT(OUT)    :: nOtherAtoms
      INTEGER,ALLOCATABLE, INTENT(OUT)    :: atomTypepList(:) !Which other atomtypes were added (not equal to refAtom)
      integer,optional,    intent(in)     :: atomicNumberSelection(:) !Which other atom should be considered

      INTEGER :: generatedShells, ishellAtom, ishell, i_gf, repr, repr_ops, atomTypep
      REAL    :: repr_diff(3)

      REAL,    ALLOCATABLE :: shellDistances(:)
      REAL,    ALLOCATABLE :: shellDiffs(:,:,:)
      INTEGER, ALLOCATABLE :: shellAtoms(:,:,:)
      INTEGER, ALLOCATABLE :: shellOps(:,:)
      INTEGER, ALLOCATABLE :: numAtomsShell(:)

      CALL timestart("Green's Function: Add nearest Neighbors")

      CALL construct_atom_shells(refAtom, nshells, atoms, cell, sym, input%film, shellDistances,&
                                 shellDiffs, shellAtoms, shellOps, numAtomsShell, generatedShells, only_elements=atomicNumberSelection)

      ALLOCATE(atomTypepList(atoms%ntype),source=0)
      nOtherAtoms = 0

      DO ishell = 1, generatedShells
         IF(l_write) THEN
            WRITE(oUnit,'(/,A,f14.8)') 'Adding shell with distance: ', shellDistances(ishell)

            WRITE(oUnit,'(/,A)') ' Contains the following atom pairs:'
            DO ishellAtom = 1, numAtomsShell(ishell)
               WRITE(oUnit,'(3f14.8,i10)') shellDiffs(:,ishellAtom,ishell), shellOps(ishellAtom,ishell)
            ENDDO
         ENDIF

         repr = 0
         DO ishellAtom = 1, numAtomsShell(ishell)

            atomTypep = atoms%itype(shellAtoms(2,ishellAtom,ishell))
            i_gf =  this%add(l,refAtom,iContour,l_sphavg,lp=lp,atomTypep=atomTypep,k_resolved=l_kresolved,&
                             atomDiff=shellDiffs(:,ishellAtom,ishell),l_fixedCutoffset=l_fixedCutoffset,&
                             fixedCutoff=fixedCutoff,atom=shellAtoms(1,ishellAtom,ishell),atomp = shellAtoms(2,ishellAtom,ishell))
            IF(repr == 0) THEN
               repr = i_gf
               repr_diff = shellDiffs(:,ishellAtom,ishell)
            ENDIF

            this%elem(i_gf)%refCutoff = refCutoff

            IF(ishellAtom > 1) THEN
               this%elem(i_gf)%representative_elem = repr
               this%elem(i_gf)%representative_op = shellops(ishellAtom,ishell)
               this%elem(i_gf)%representative_diff = repr_diff
            ENDIF

            IF(atomTypep.NE.refAtom.AND..NOT.ANY(atomTypepList(:nOtherAtoms).EQ.atomTypep)) THEN
               !Other atomtype
               nOtherAtoms = nOtherAtoms + 1
               atomTypepList(nOtherAtoms) = atomTypep
            ENDIF


            IF(l_write) THEN
               WRITE(oUnit,'(A,I6,I6,3f14.8,3i10)') 'GF Element: ', refAtom, atomTypep,&
                                               shellDiffs(:,ishellAtom,ishell), &
                                               shellops(ishellAtom,ishell), &
                                               shellAtoms(:,ishellAtom,ishell)
            ENDIF
         ENDDO

      ENDDO


   CALL timestop("Green's Function: Add nearest Neighbors")

   END SUBROUTINE addNearestNeighbours_gfelem

   INTEGER FUNCTION find_gfelem_simple(this,l,atomType,iContour,l_sphavg,lp,atomTypep,&
                                       nshells,atomDiff,k_resolved,k_resolved_int,l_found) result(i_gf)

      !Maps between the four indices (l,lp,nType,nTypep) and the position in the
      !gf arrays

      CLASS(t_gfinp),      INTENT(IN)    :: this
      INTEGER,             INTENT(IN)    :: l
      INTEGER,             INTENT(IN)    :: atomType
      INTEGER,             INTENT(IN)    :: iContour
      LOGICAL,             INTENT(IN)    :: l_sphavg
      INTEGER, OPTIONAL,   INTENT(IN)    :: lp
      INTEGER, OPTIONAL,   INTENT(IN)    :: atomTypep
      INTEGER, OPTIONAL,   INTENT(IN)    :: nshells
      REAL,    OPTIONAL,   INTENT(IN)    :: atomDiff(:)
      LOGICAL, OPTIONAL,   INTENT(IN)    :: k_resolved
      LOGICAL, OPTIONAL,   INTENT(IN)    :: k_resolved_int
      LOGICAL, OPTIONAL,   INTENT(INOUT) :: l_found    !If this switch is not provided the program
                                                       !will assume that the element has to be present and
                                                       !terminate with an error message if the
                                                       !element is not found (for adding elements)

      LOGICAL :: search, distinct_kresolved_int
      TYPE(t_gfelementtype) :: elem_to_find

      distinct_kresolved_int = PRESENT(k_resolved_int)

      CALL elem_to_find%init(l,atomType,iContour,l_sphavg,lp=lp,atomTypep=atomTypep,&
                             nshells=nshells,atomDiff=atomDiff,k_resolved=k_resolved)
      IF(distinct_kresolved_int) THEN
         elem_to_find%l_kresolved_int = k_resolved_int
      ENDIF

      i_gf = this%find(elem_to_find, l_found=l_found, &
                       distinct_kresolved_int=distinct_kresolved_int)

   END FUNCTION find_gfelem_simple

   INTEGER FUNCTION find_gfelem_type(this,elem,l_found, distinct_kresolved_int) result(i_gf)

      !Maps between the four indices (l,lp,nType,nTypep) and the position in the
      !gf arrays

      CLASS(t_gfinp),         INTENT(IN)    :: this
      TYPE(t_gfelementtype),  INTENT(IN)    :: elem
      LOGICAL, OPTIONAL,      INTENT(INOUT) :: l_found    !If this switch is not provided the program
                                                          !will assume that the element has to be present and
                                                          !terminate with an error message if the
                                                          !element is not found (for adding elements)
      LOGICAL, OPTIONAL,      INTENT(IN)    :: distinct_kresolved_int

      LOGICAL :: search

      IF(PRESENT(l_found)) l_found = .FALSE.
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
         IF(.NOT.elem%equals(this%elem(i_gf), distinct_k_resolved=distinct_kresolved_int)) CYCLE

         !If we are here we found the element
         IF(PRESENT(l_found)) l_found=.TRUE.
         search = .FALSE.
      ENDDO
   END FUNCTION find_gfelem_type

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
         IF(.NOT.this%elem(i_gf)%isOffDiag()) checkOnsite_gfinp = .TRUE.
      ENDDO
   END FUNCTION checkOnsite_gfinp

   PURE LOGICAL FUNCTION checkOffdiagonal_gfinp(this)

      !Check if there are any oniste elements
      CLASS(t_gfinp),               INTENT(IN)    :: this

      INTEGER :: i_gf

      checkOffdiagonal_gfinp = .FALSE.
      DO i_gf = 1, this%n
         IF(this%elem(i_gf)%isOffDiag()) checkOffdiagonal_gfinp = .TRUE.
      ENDDO
   END FUNCTION checkOffdiagonal_gfinp

   SUBROUTINE init_gfelem(this,l,atomType,iContour,l_sphavg,lp,nshells,atomTypep,k_resolved,atomDiff,l_fixedCutoffset,fixedCutoff,atom,atomp)

      CLASS(t_gfelementtype), INTENT(INOUT)  :: this
      INTEGER,                INTENT(IN)     :: l
      INTEGER,                INTENT(IN)     :: atomType
      INTEGER,                INTENT(IN)     :: iContour
      LOGICAL,                INTENT(IN)     :: l_sphavg
      INTEGER, OPTIONAL,      INTENT(IN)     :: lp
      INTEGER, OPTIONAL,      INTENT(IN)     :: atomTypep !Specify the second atom
      REAL,    OPTIONAL,      INTENT(IN)     :: atomDiff(:)
      LOGICAL, OPTIONAL,      INTENT(IN)     :: l_fixedCutoffset
      REAL,    OPTIONAL,      INTENT(IN)     :: fixedCutoff
      INTEGER, OPTIONAL,      INTENT(IN)     :: nshells
      LOGICAL, OPTIONAL,      INTENT(IN)     :: k_resolved
      INTEGER, OPTIONAL,      INTENT(IN)     :: atom, atomp

      IF(PRESENT(nshells).AND.PRESENT(atomTypep)) THEN
         CALL juDFT_error("Conflicting arguments: nshells and nTypep given",&
                          hint="This is a bug in FLEUR, please report",&
                          calledby="init_gfelem")
      ENDIF

      IF(PRESENT(atom).OR.PRESENT(atomp)) THEN
         IF(.NOT.(PRESENT(atom).AND.PRESENT(atomp))) THEN
            CALL juDFT_error("Invalid arguments: Either both atom and atomp need to be given or none of them",&
                             hint="This is a bug in FLEUR, please report",&
                             calledby="init_gfelem")
         ENDIF
      ENDIF

      this%l = l
      this%atomType = atomType
      this%iContour = iContour
      this%l_sphavg = l_sphavg
      IF(PRESENT(lp)) THEN
         this%lp = lp
      ELSE
         this%lp = l
      ENDIF
      IF(PRESENT(nshells)) THEN
         IF(nshells/=0) THEN
            !Temporary index to mark later in gfinp%init
            this%atomTypep = -nshells
         ELSE
            this%atomTypep = atomType
         ENDIF
      ELSE IF(PRESENT(atomTypep)) THEN
         !Explicit declaration of element
         this%atomTypep = atomTypep
      ELSE
         !No intersite element
         this%atomTypep = atomType
      ENDIF
      IF(PRESENT(k_resolved)) THEN
         this%l_kresolved = k_resolved
         this%l_kresolved_int = k_resolved
      ENDIF
      IF(PRESENT(atomDiff)) THEN
         this%atomDiff(:) = atomDiff(:)
      ELSE
         this%atomDiff(:) = 0.0
      ENDIF

      IF(PRESENT(atom)) THEN
         this%atom = atom
         this%atomp = atomp
      ELSE
         this%atom = 0
         this%atomp = 0
      ENDIF

      IF(PRESENT(l_fixedCutoffset)) THEN
         IF(.NOT.PRESENT(fixedCutoff)) CALL juDFT_error("l_fixedCutoffset Present without fixedCutoff", &
                                                        hint="This is a bug in FLEUR please report",&
                                                        calledby="init_gfelem")
         this%l_fixedCutoffset = l_fixedCutoffset
         IF(l_fixedCutoffset) THEN
            this%fixedCutoff = fixedCutoff
         ENDIF
      ENDIF

   END SUBROUTINE init_gfelem

   PURE LOGICAL FUNCTION equals_coefficients_gfelem(this, other, distinct_k_resolved, distinct_symmetry_equivalent_diffs)

      CLASS(t_gfelementtype), INTENT(IN)  :: this
      TYPE(t_gfelementtype),  INTENT(IN)  :: other
      LOGICAL, OPTIONAL,      INTENT(IN)  :: distinct_k_resolved
      LOGICAL, OPTIONAL,      INTENT(IN)  :: distinct_symmetry_equivalent_diffs

      LOGICAL distinct_k_resolved_arg, distinct_symmetry_equivalent_diffs_arg

      distinct_k_resolved_arg = .TRUE.
      IF(PRESENT(distinct_k_resolved)) distinct_k_resolved_arg = distinct_k_resolved

      distinct_symmetry_equivalent_diffs_arg = .FALSE.
      IF(PRESENT(distinct_symmetry_equivalent_diffs)) distinct_symmetry_equivalent_diffs_arg = distinct_symmetry_equivalent_diffs

      equals_coefficients_gfelem = .FALSE.

      IF(this%l.NE.other%l) RETURN
      IF(this%lp.NE.other%lp) RETURN
      IF(this%atomType.NE.other%atomType) RETURN
      IF(this%atomTypep.NE.other%atomTypep) RETURN
      IF(this%l_sphavg .neqv. other%l_sphavg) RETURN
      IF(this%l_kresolved .neqv. other%l_kresolved) RETURN
      IF(distinct_k_resolved_arg) then
         IF(this%l_kresolved_int .neqv. other%l_kresolved_int) RETURN
      ENDIF
      IF(ANY(ABS(this%atomDiff(:)-other%atomDiff(:)).GT.ATOMDIFF_EPS)) THEN
         IF(distinct_symmetry_equivalent_diffs_arg) RETURN
         IF(this%representative_elem < 0 .AND. other%representative_elem < 0) RETURN
         IF(this%representative_elem > 0 .AND. other%representative_elem > 0) THEN
            IF(ANY(ABS(this%representative_diff(:)-other%representative_diff(:)).GT.ATOMDIFF_EPS)) RETURN
         ELSE IF(this%representative_elem > 0) THEN
            IF(ANY(ABS(this%representative_diff(:)-other%atomDiff(:)).GT.ATOMDIFF_EPS)) RETURN
         ELSE IF(other%representative_elem > 0) THEN
            IF(ANY(ABS(other%representative_diff(:)-this%atomDiff(:)).GT.ATOMDIFF_EPS)) RETURN
         ENDIF
      ENDIF
      IF(this%atom /= other%atom) RETURN
      IF(this%atomp /= other%atomp) RETURN
      equals_coefficients_gfelem = .TRUE.

   END FUNCTION equals_coefficients_gfelem

   PURE LOGICAL FUNCTION equals_gfelem(this, other, distinct_k_resolved)

      CLASS(t_gfelementtype), INTENT(IN)  :: this
      TYPE(t_gfelementtype),  INTENT(IN)  :: other
      LOGICAL, OPTIONAL,      INTENT(IN)  :: distinct_k_resolved

      equals_gfelem = .FALSE.
      !We need to check the atomDiff explicitly, since the deduplication
      !on the coefficient level has some extra symmetry considerations
      !that should not influence the deduplication. It just influences how
      !many brillouin zone integegrations need to be performed
      IF(.NOT.this%equals_coefficients(other, distinct_k_resolved, distinct_symmetry_equivalent_diffs=.TRUE.)) RETURN
      IF(this%iContour.NE.other%iContour) RETURN
      equals_gfelem = .TRUE.

   END FUNCTION equals_gfelem

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

      isoffDiag_gfelem = this%l.NE.this%lp.OR.this%isIntersite()

   END FUNCTION isOffDiag_gfelem

   PURE LOGICAL FUNCTION isIntersite_gfelem(this)

      CLASS(t_gfelementtype),   INTENT(IN)  :: this

      isIntersite_gfelem = this%atomType.NE.this%atomTypep&
                           .OR.ANY(ABS(this%atomDiff).GT.ATOMDIFF_EPS)

   END FUNCTION isIntersite_gfelem

END MODULE m_types_gfinp
