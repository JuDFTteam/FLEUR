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
      LOGICAL :: l_fixedCutoffset = .FALSE.
      REAL    :: fixedCutoff = 0.0
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
      LOGICAL :: l_sphavg = .TRUE.
      LOGICAL :: l_mperp = .FALSE.
      !Number of elements !(TODO: NAME??)
      INTEGER :: n = 0
      !Number of j0 calculations
      INTEGER :: n_j0 = 0
      !Information on the elements to be calculated
      TYPE(t_gfelementtype), ALLOCATABLE :: elem(:)
      !Parameters for the energy mesh on the real axis
      INTEGER :: ne    = 2700
      REAL    :: ellow = -1.0
      REAL    :: elup  =  1.0
      INTEGER :: numberContours = 0
      TYPE(t_contourInp), ALLOCATABLE :: contour(:)
      INTEGER, ALLOCATABLE :: hiaElem(:)
   CONTAINS
      PROCEDURE :: read_xml      => read_xml_gfinp
      PROCEDURE :: mpi_bc        => mpi_bc_gfinp
      PROCEDURE :: init          => init_gfinp
      PROCEDURE :: find          => find_gfelem
      PROCEDURE :: find_contour  => find_contour
      PROCEDURE :: add           => add_gfelem
      PROCEDURE :: eMesh         => eMesh_gfinp
      PROCEDURE :: addNearestNeighbours => addNearestNeighbours_gfelem
   END TYPE t_gfinp

   PUBLIC t_gfinp, t_contourInp
   PUBLIC CONTOUR_RECTANGLE_CONST, CONTOUR_SEMICIRCLE_CONST, CONTOUR_DOS_CONST
   PUBLIC uniqueElements_gfinp

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
      CALL mpi_bc(this%l_sphavg,rank,mpi_comm)
      CALL mpi_bc(this%l_mperp,rank,mpi_comm)
      CALL mpi_bc(this%n,rank,mpi_comm)
      CALL mpi_bc(this%n_j0,rank,mpi_comm)
      CALL mpi_bc(this%ne,rank,mpi_comm)
      CALL mpi_bc(this%ellow,rank,mpi_comm)
      CALL mpi_bc(this%elup,rank,mpi_comm)
      CALL mpi_bc(this%numberContours,rank,mpi_comm)
      CALL mpi_bc(this%hiaElem,rank,mpi_comm)

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

      INTEGER :: numberNodes,ntype,itype,n_hia
      INTEGER :: lmin,lmax,i,l,lp,iContour,iContourp
      REAL    :: fixedCutoff
      CHARACTER(len=100)  :: xPathA,xPathS,label,cutoffArg
      LOGICAL :: l_gfinfo_given,l_off,l_nn,l_fixedCutoffset

      xPathA = '/fleurInput/calculationSetup/greensFunction'
      numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
      l_gfinfo_given = numberNodes.EQ.1

      IF (l_gfinfo_given) THEN
         this%l_sphavg=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_sphavg'))
         this%l_mperp=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))

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
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/contourRectangle[',i,']'
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
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/contourSemicircle[',i,']'
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
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/contourDOS[',i,']'
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

      ALLOCATE(this%elem(4*ntype))
      ALLOCATE(this%hiaElem(4*ntype))

      DO itype = 1, ntype
         xPathS=xml%speciesPath(itype)

         !Read in all possible tags, which need a greens function calculation

         !Explicit declaration of an onsite greens function element
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/onsiteGF')
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/onsiteGF[',i,']'
            lmin = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_min'))
            lmax = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_max'))
            l_off = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_offdiag'))
            label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
            iContour = this%find_contour(TRIM(ADJUSTL(label)))
            cutoffArg = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kkintgrCutoff')))
            IF(TRIM(ADJUSTL(cutoffArg))=="calc") THEN
               l_fixedCutoffset = .FALSE.
            ELSE
               fixedCutoff = evaluateFirstOnly(TRIM(ADJUSTL(cutoffArg)))
            ENDIF
            DO l = lmin, lmax
               DO lp = MERGE(lmin,l,l_off), MERGE(lmax,l,l_off)
                  CALL this%add(itype,l,lp,iContour,l_fixedCutoffset=l_fixedCutoffset,&
                                fixedCutoff=fixedCutoff)
               ENDDO
            ENDDO
         ENDDO

         !Explicit declaration of an intersite greens function element
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/intersiteGF')
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/intersiteGF[',i,']'
            lmin = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_min'))
            lmax = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_max'))
            l_off = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_offdiag'))
            label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
            cutoffArg = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kkintgrCutoff')))
            IF(TRIM(ADJUSTL(cutoffArg))=="calc") THEN
               l_fixedCutoffset = .FALSE.
            ELSE
               fixedCutoff = evaluateFirstOnly(TRIM(ADJUSTL(cutoffArg)))
            ENDIF
            iContour = this%find_contour(TRIM(ADJUSTL(label)))
            DO l = lmin, lmax
               DO lp = MERGE(lmin,l,l_off), MERGE(lmax,l,l_off)
                  CALL this%add(itype,l,lp,iContour,l_fixedCutoffset=l_fixedCutoffset,&
                                fixedCutoff=fixedCutoff,l_inter=.TRUE.)
               ENDDO
            ENDDO
         ENDDO

         !Declaration of a DFT+Hubbard 1 calculation
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/ldaHIA')
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/ldaHIA[',i,']'
            !No offdiagonal l-part (only a single element)
            l = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l'))
            label = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@label')))
            iContour = this%find_contour(TRIM(ADJUSTL(label)))
            cutoffArg = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kkintgrCutoff')))
            IF(TRIM(ADJUSTL(cutoffArg))=="calc") THEN
               l_fixedCutoffset = .FALSE.
            ELSE
               fixedCutoff = evaluateFirstOnly(TRIM(ADJUSTL(cutoffArg)))
            ENDIF
            CALL this%add(itype,l,l,iContour,l_fixedCutoffset=l_fixedCutoffset,&
                          fixedCutoff=fixedCutoff)
            n_hia = n_hia + 1
            this%hiaElem(n_hia) = this%n
         ENDDO
      ENDDO

      IF(this%n>0.AND..NOT.l_gfinfo_given) THEN
         CALL juDFT_error("Error reading in gf-information: No general information found for the gf-calculations",&
                           calledby="read_xml_gfinp")
      ENDIF

      !Check the input for validity
      IF(this%n.GT.0) THEN
         IF(this%elup.GT.1.0) CALL juDFT_warn("Cutoff for the Greens function calculation should never be higher"//&
                                              "than 1htr above efermi",calledby="read_xml_gfinp")
         IF(this%elup.LT.this%ellow) CALL juDFT_error("Not a valid energy grid elup<ellow",calledby="read_xml_gfinp")
         IF(ANY(this%elem(:this%n)%l.LT.2)) CALL juDFT_warn("Green's function for s and p orbitals not tested",&
                                                            calledby="read_xml_gfinp")
         IF(ANY(this%elem(:this%n)%l.GT.3)) CALL juDFT_error("Green's function only implemented for l<=3",&
                                                             calledby="read_xml_gfinp")
      ENDIF

   END SUBROUTINE read_xml_gfinp

   SUBROUTINE init_gfinp(this,atoms,sym,noco)

      USE m_types_atoms
      USE m_types_sym
      USE m_types_noco

      CLASS(t_gfinp),   INTENT(INOUT)  :: this
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_noco),     INTENT(IN)     :: noco

      INTEGER :: i_gf,l,lp,atomType,atomTypep,iContour
      INTEGER :: hiaElem(atoms%n_hia)

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
                                           this%elem(i_gf)%fixedCutoff,atoms,sym)
         ENDIF
      ENDDO

      !Reallocate with correct size
      hiaElem = this%hiaElem(:atoms%n_hia)
      IF(ALLOCATED(this%hiaElem)) DEALLOCATE(this%hiaElem)
      ALLOCATE(this%hiaElem(atoms%n_hia))
      this%hiaElem = hiaElem

      IF(this%l_mperp.AND..NOT.noco%l_mperp) THEN
         CALL juDFT_error("For l_mperp for Green's Functions the l_mperp switch for noco has to be True",&
                          calledby="init_gfinp")
      ENDIF

   END SUBROUTINE init_gfinp

   SUBROUTINE uniqueElements_gfinp(gfinp,uniqueElements,ind,indUnique)

      !Not a procedure, because gfortran+OpenMP has problems with it
      !Called inside OMP parallel region

      TYPE(t_gfinp),    INTENT(IN)     :: gfinp
      INTEGER,          INTENT(INOUT)  :: uniqueElements !Number of unique elements before ind or in the whole array
      INTEGER, OPTIONAL,INTENT(IN)     :: ind
      INTEGER, OPTIONAL,INTENT(INOUT)  :: indUnique      !Position of the corresponding unique Element for a given ind

      INTEGER :: maxGF
      INTEGER :: l,lp,atomType,atomTypep,dummyInd,iContour,i_gf
      LOGICAL :: l_unique

      uniqueElements = 0

      IF(PRESENT(ind)) THEN
         maxGF = ind
      ELSE
         maxGF = gfinp%n
      ENDIF
      DO i_gf = 1, maxGF
         l  = gfinp%elem(i_gf)%l
         lp = gfinp%elem(i_gf)%lp
         atomType  = gfinp%elem(i_gf)%atomType
         atomTypep = gfinp%elem(i_gf)%atomTypep
         iContour  = gfinp%elem(i_gf)%iContour
         dummyInd = gfinp%find(l,atomType,iContour=iContour,lp=lp,nTypep=atomTypep,&
                              uniqueMax=i_gf,l_unique=l_unique)
         IF(l_unique) THEN
            uniqueElements = uniqueElements +1
         ENDIF
      ENDDO

      IF(uniqueElements==0.AND.maxGF/=0) THEN
         CALL juDFT_error("No unique GF elements",hint="This is a bug in FLEUR please report",&
                          calledby="uniqueElements_gfinp")
      ENDIF

      IF(PRESENT(indUnique)) THEN
         IF(.NOT.PRESENT(ind)) CALL juDFT_error("ind and indUnique have to be provided at the same time",&
                                                calledby="uniqueElements_gfinp")
         l  = gfinp%elem(ind)%l
         lp = gfinp%elem(ind)%lp
         atomType  = gfinp%elem(ind)%atomType
         atomTypep = gfinp%elem(ind)%atomTypep
         iContour  = gfinp%elem(ind)%iContour

         indUnique = gfinp%find(l,atomType,iContour=iContour,lp=lp,nTypep=atomTypep,&
                               uniqueMax=ind,l_unique=l_unique)
      ENDIF

   END SUBROUTINE uniqueElements_gfinp

   SUBROUTINE add_gfelem(this,nType,l,lp,iContour,nTypep,l_fixedCutoffset,fixedCutoff,l_inter)

      CLASS(t_gfinp),      INTENT(INOUT)  :: this
      INTEGER,             INTENT(IN)     :: nType
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: lp
      INTEGER,             INTENT(IN)     :: iContour
      INTEGER, OPTIONAL,   INTENT(IN)     :: nTypep !Specify the second atom
      LOGICAL, OPTIONAL,   INTENT(IN)     :: l_fixedCutoffset
      REAL,    OPTIONAL,   INTENT(IN)     :: fixedCutoff
      LOGICAL, OPTIONAL,   INTENT(IN)     :: l_inter!To be used in init when atoms is not available and nTypep was not specified


      INTEGER i_gf
      LOGICAL l_found

      IF(PRESENT(l_inter).AND.PRESENT(nTypep)) CALL juDFT_error("Conflicting arguments: l_inter and nTypep given",&
                                                                hint="This is a bug in FLEUR, please report",&
                                                                calledby="add_gfelem")

      !Check if this job has already been added
      i_gf = this%find(l,nType,lp,nType,iContour=iContour,l_found=l_found)
      IF(l_found) RETURN !Element was found

      this%n = this%n + 1
      this%elem(this%n)%l = l
      this%elem(this%n)%atomType = nType
      this%elem(this%n)%lp = lp
      this%elem(this%n)%iContour = iContour
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

   END SUBROUTINE add_gfelem

   SUBROUTINE addNearestNeighbours_gfelem(this,nshells,l,lp,refAtom,iContour,l_fixedCutoffset,fixedCutoff,atoms,sym)

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
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym

      INTEGER :: ishell,natomp
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
               CALL this%add(refAtom,l,lp,iContour,nTypep=natomp,l_fixedCutoffset=l_fixedCutoffset,&
                             fixedCutoff=fixedCutoff)
               dist(natomp) = 9e99 !Eliminate from the list
            ENDIF
         ENDDO
         ishell = ishell + 1
      ENDDO

   END SUBROUTINE addNearestNeighbours_gfelem

   FUNCTION find_gfelem(this,l,nType,lp,nTypep,iContour,uniqueMax,l_unique,l_found) result(i_gf)

      !Maps between the four indices (l,lp,nType,nTypep) and the position in the
      !gf arrays

      CLASS(t_gfinp),      INTENT(IN)    :: this
      INTEGER,             INTENT(IN)    :: l
      INTEGER,             INTENT(IN)    :: nType
      INTEGER, OPTIONAL,   INTENT(IN)    :: iContour
      INTEGER, OPTIONAL,   INTENT(IN)    :: lp
      INTEGER, OPTIONAL,   INTENT(IN)    :: nTypep

      INTEGER, OPTIONAL,   INTENT(IN)    :: uniqueMax  !These arguments will return whether there
      LOGICAL, OPTIONAL,   INTENT(INOUT) :: l_unique   !is an element before uniqueMax with the same (l,lp,nType,nTypep)
                                                       !combination but different energy contour

      LOGICAL, OPTIONAL,   INTENT(INOUT) :: l_found    !If this switch is not provided the program
                                                       !will assume that the element has to be present and
                                                       !terminate with an error message if the
                                                       !element is not found (for adding elements)

      INTEGER :: i_gf
      LOGICAL :: search

      search = .TRUE.
      IF(PRESENT(l_unique)) l_unique = .TRUE.
      IF((PRESENT(l_unique).OR.PRESENT(uniqueMax)).AND..NOT.PRESENT(l_unique).AND.PRESENT(uniqueMax)) THEN
         CALL juDFT_error("Not provided uniqueMax AND l_unique",&
                          hint="This is a bug in FLEUR please report",calledby="find_gfelem")
      ENDIF
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
         !If we are here and smaller than uniqueMax the element is not unique
         IF(PRESENT(uniqueMax)) THEN
            IF(i_gf<uniqueMax) THEN
               l_unique = .FALSE.
               RETURN !Return the unique index for the element
            ENDIF
         ENDIF
         IF(PRESENT(iContour)) THEN
            IF(this%elem(i_gf)%iContour.NE.iContour) CYCLE
         ELSE
            IF(this%elem(i_gf)%iContour.NE.1) CYCLE
         ENDIF
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

END MODULE m_types_gfinp
