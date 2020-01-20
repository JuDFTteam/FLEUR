!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_gfinp
   USE m_juDFT
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE

   TYPE t_gfelementtype
      SEQUENCE
      !defines the l and atomType elements for given greens function element
      !(used for mapping index in types_greensf)
      INTEGER :: l = -1
      INTEGER :: lp = -1
      INTEGER :: atomType = 0
      INTEGER :: atomTypep = 0
   END TYPE t_gfelementtype

   TYPE t_j0calctype
      SEQUENCE
      INTEGER :: atomType = 0             !atom Type for which to calculate J0
      INTEGER :: lmin = -1                !Minimum l considered
      INTEGER :: lmax = -1                !Maximum l considered
      LOGICAL :: l_avgexc = .FALSE.       !Determines whether to average over the exchange splittings for all l
      LOGICAL :: l_eDependence = .FALSE.  !Switch to output J0 with variing fermi energy (only with contourDOS)
   END TYPE t_j0calctype


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
      !Information on the j0-elements to be calculated
      TYPE(t_j0calctype), ALLOCATABLE    :: j0elem(:)
      !Parameters for the energy mesh on the real axis
      INTEGER :: ne = 1301
      REAL    :: ellow = -1.0
      REAL    :: elup = 1.0
      INTEGER :: mode = -1 !Controls the shape of the complex energy contour
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
   CONTAINS
      PROCEDURE :: read_xml   => read_xml_gfinp
      PROCEDURE :: mpi_bc     => mpi_bc_gfinp
      PROCEDURE :: find       => find_gfelem
      PROCEDURE :: add        => add_gfelem
      PROCEDURE :: eMesh      => eMesh_gfinp
      PROCEDURE :: eContour   => eContour_gfinp
   END TYPE t_gfinp
   PUBLIC t_gfinp

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
      CALL mpi_bc(this%mode,rank,mpi_comm)
      CALL mpi_bc(this%eb,rank,mpi_comm)
      CALL mpi_bc(this%et,rank,mpi_comm)
      CALL mpi_bc(this%n1,rank,mpi_comm)
      CALL mpi_bc(this%n2,rank,mpi_comm)
      CALL mpi_bc(this%n3,rank,mpi_comm)
      CALL mpi_bc(this%nmatsub,rank,mpi_comm)
      CALL mpi_bc(this%sigma,rank,mpi_comm)
      CALL mpi_bc(this%ncirc,rank,mpi_comm)
      CALL mpi_bc(this%alpha,rank,mpi_comm)
      CALL mpi_bc(this%nDOS,rank,mpi_comm)
      CALL mpi_bc(this%sigmaDOS,rank,mpi_comm)
      CALL mpi_bc(this%l_anacont,rank,mpi_comm)
      CALL mpi_bc(this%l_dosfermi,rank,mpi_comm)

#ifdef CPP_MPI
      CALL mpi_COMM_RANK(mpi_comm,myrank,ierr)
      IF (myrank.NE.rank) THEN
         IF (ALLOCATED(this%elem)) DEALLOCATE(this%elem)
         IF (ALLOCATED(this%j0elem)) DEALLOCATE(this%j0elem)
         ALLOCATE(this%elem(this%n))
         ALLOCATE(this%j0elem(this%n_j0))
      ENDIF
      DO n=1,this%n
         CALL mpi_bc(this%elem(n)%l,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%atomType,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%lp,rank,mpi_comm)
         CALL mpi_bc(this%elem(n)%atomTypep,rank,mpi_comm)
      ENDDO
      DO n=1,this%n_j0
         CALL mpi_bc(this%j0elem(n)%atomType,rank,mpi_comm)
         CALL mpi_bc(this%j0elem(n)%lmin,rank,mpi_comm)
         CALL mpi_bc(this%j0elem(n)%lmax,rank,mpi_comm)
         CALL mpi_bc(this%j0elem(n)%l_avgexc,rank,mpi_comm)
         CALL mpi_bc(this%j0elem(n)%l_eDependence,rank,mpi_comm)
      ENDDO
#endif

   END SUBROUTINE mpi_bc_gfinp

   SUBROUTINE read_xml_gfinp(this, xml)
      USE m_types_xml
      CLASS(t_gfinp), INTENT(INOUT):: this
      TYPE(t_xml),INTENT(INOUT) ::xml

      INTEGER :: numberNodes,ntype,itype
      INTEGER :: lmin,lmax,i
      CHARACTER(len=100)  :: xPathA,xPathS
      LOGICAL :: l_gfinfo_given,l_off,l_nn

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

         !Information on the complex energy contour
         xPathA = '/fleurInput/calculationSetup/greensFunction/contourRectangle'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF(numberNodes.EQ.1) THEN
            this%mode = 1
            this%eb = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eb'))
            !et cannot be varied from the fermi energy for this contour
            this%n1 = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n1'))
            this%n2 = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n2'))
            this%n3 = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n3'))
            this%nmatsub = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nmatsub'))
            this%sigma = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sigma'))
         ENDIF

         xPathA = '/fleurInput/calculationSetup/greensFunction/contourSemicircle'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF(numberNodes.EQ.1) THEN
            this%mode = 2
            this%eb = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eb'))
            this%et = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@et'))
            this%ncirc = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n'))
            this%alpha = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@alpha'))
         ENDIF

         xPathA = '/fleurInput/calculationSetup/greensFunction/contourDOS'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF(numberNodes.EQ.1) THEN
            this%mode = 3
            this%eb = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eb'))
            this%et = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@et'))
            this%nDOS = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@n'))
            this%sigmaDOS = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sigma'))
            this%l_anacont = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@analytical_cont'))
            this%l_dosfermi = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_fermi'))
         ENDIF

         IF(this%mode.EQ.-1) CALL juDFT_error("Error reading in gf-information: No complex energy contour specified",&
                                              calledby="read_xml_gfinp")
      ENDIF

      ntype = xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')

      ALLOCATE(this%elem(4*ntype))
      ALLOCATE(this%j0elem(4*ntype))

      DO itype = 1, ntype
         xPathS=xml%speciesPath(itype)

         !Read in all possible tags, which need a greens function calculation

         !Explicit declaration of an onsite greens function element
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/onsiteGF')
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/onsiteGF[',i,']'
            lmin = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_min'))
            lmax = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_max'))
            l_off = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_offdiag'))
            CALL this%add(itype,lmin,lmax,l_off,.FALSE.,.FALSE.)
         ENDDO

         !Explicit declaration of an intersite greens function element
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/intersiteGF')
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/intersiteGF[',i,']'
            lmin = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_min'))
            lmax = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_max'))
            l_off = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_offdiag'))
            l_nn = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_nn'))
            CALL this%add(itype,lmin,lmax,l_off,.TRUE.,l_nn)
         ENDDO

         !Declaration of a j0 calculation
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/J0')
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/J0[',i,']'
            lmin = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_min'))
            lmax = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_max'))
            CALL this%add(itype,lmin,lmax,.FALSE.,.FALSE.,.FALSE.)
            !Add it to the j0elem array
            this%n_j0 = this%n_j0 + 1
            this%j0elem(this%n_j0)%atomType = itype
            this%j0elem(this%n_j0)%lmin     = lmin
            this%j0elem(this%n_j0)%lmax     = lmax
            this%j0elem(this%n_j0)%l_avgexc = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_avgexc'))
            this%j0elem(this%n_j0)%l_eDependence = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_eDependence'))
         ENDDO

         !Declaration of a DFT+Hubbard 1 calculation
         DO i = 1, xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/ldaHIA')
            WRITE(xPathA,*) TRIM(ADJUSTL(xPathS))//'/ldaHIA[',i,']'
            !No offdiagonal l-part (only a single element)
            lmin = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l'))
            CALL this%add(itype,lmin,lmin,.FALSE.,.FALSE.,.FALSE.)
         ENDDO
      ENDDO

      IF(this%n>0.AND..NOT.l_gfinfo_given) CALL juDFT_error("Error reading in gf-information: No general information found for the gf-calculations",&
                                              calledby="read_xml_gfinp")

      !Check the input for validity
      IF(this%n.GT.0) THEN
         IF(this%elup.GT.1.0) CALL juDFT_warn("Cutoff for the Greens function calculation should never be higher"//&
                                              "than 1htr above efermi",calledby="read_xml_gfinp")
         IF(this%elup.LT.this%ellow) CALL juDFT_error("Not a valid energy grid elup<ellow",calledby="read_xml_gfinp")
         IF(ANY(this%elem(:this%n)%l.LT.2)) CALL juDFT_warn("Green's function for s and p orbitals not tested",calledby="read_xml_gfinp")
         IF(ANY(this%elem(:this%n)%l.GT.3)) CALL juDFT_error("Green's function only implemented for l<=3",calledby="read_xml_gfinp")

         DO i = 1, this%n_j0
            IF(this%j0elem(i)%lmin.GT.this%j0elem(i)%lmax) CALL juDFT_error("Not a valid configuration for J0-calculation l_min>l_max", &
                                                                    calledby="read_xml_gfinp")
            IF(this%j0elem(i)%l_eDependence.AND.this%mode.NE.3) CALL juDFT_error("Energy dependence of J0 only available with contourDOS",&
                                                                            calledby="read_xml_gfinp")
         ENDDO
      ENDIF

   END SUBROUTINE read_xml_gfinp

   SUBROUTINE add_gfelem(this,nType,lmin,lmax,l_off,l_inter,l_nn)

      CLASS(t_gfinp),   INTENT(INOUT)  :: this
      INTEGER,          INTENT(IN)     :: nType
      INTEGER,          INTENT(IN)     :: lmin
      INTEGER,          INTENT(IN)     :: lmax
      LOGICAL,          INTENT(IN)     :: l_off !l!=lp
      LOGICAL,          INTENT(IN)     :: l_inter
      LOGICAL,          INTENT(IN)     :: l_nn

      INTEGER l,lp,i_gf
      LOGICAL l_found

      IF(l_inter) CALL juDFT_error("Intersite greens function not yet implemented",calledby="add_gfjob")

      !TODO: add the nearest neighbours jobs

      DO l = lmin, lmax
         DO lp = MERGE(lmin,l,l_off), MERGE(lmax,l,l_off)
            !Check if this job has already been added (notice 2 ntype) To be added!!
            i_gf = this%find(l,nType,lp,nType,l_found)
            IF(l_found) CYCLE !Element was found

            this%n = this%n + 1
            this%elem(this%n)%l = l
            this%elem(this%n)%atomType = nType
            this%elem(this%n)%lp = lp
            this%elem(this%n)%atomTypep = nType !For now

         ENDDO
      ENDDO

   END SUBROUTINE add_gfelem

   FUNCTION find_gfelem(this,l,nType,lp,nTypep,l_found) result(i_gf)

      !Maps between the four indices (l,lp,nType,nTypep) and the position in the
      !gf arrays

      CLASS(t_gfinp),      INTENT(IN)    :: this
      INTEGER,             INTENT(IN)    :: l
      INTEGER,             INTENT(IN)    :: nType
      INTEGER, OPTIONAL,   INTENT(IN)    :: lp
      INTEGER, OPTIONAL,   INTENT(IN)    :: nTypep
      LOGICAL, OPTIONAL,   INTENT(INOUT) :: l_found !If this switch is not provided the program
                                                  !will terminate with an error message if the
                                                  !element is not found (for adding elements)

      INTEGER :: i_gf
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
            IF(this%elem(i_gf)%atomTypep.NE.nTypep) CYCLE
         ELSE
            IF(this%elem(i_gf)%atomTypep.NE.nType) CYCLE
         ENDIF
         !If we are here we found the element
         IF(PRESENT(l_found)) l_found=.TRUE.
         search = .FALSE.
      ENDDO
   END FUNCTION find_gfelem

   SUBROUTINE eMesh_gfinp(this,ef,del_out,eb_out,et_out,eMesh)

      !Gives back the information for the energy mesh on the real axis
      !Energies are shifted according to the fermi level

      CLASS(t_gfinp),               INTENT(IN)    :: this
      REAL,                         INTENT(IN)    :: ef
      REAL, OPTIONAL,               INTENT(INOUT) :: del_out
      REAL, OPTIONAL,               INTENT(INOUT) :: eb_out
      REAL, OPTIONAL,               INTENT(INOUT) :: et_out
      REAL, ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: eMesh(:)

      INTEGER :: ie
      REAL :: del,eb

      eb = ef+this%ellow
      del = (this%elup-this%ellow)/REAL(this%ne-1)

      IF(PRESENT(eb_out)) eb_out = eb
      IF(PRESENT(et_out)) et_out = ef+this%elup
      IF(PRESENT(del_out)) del_out = del

      IF(PRESENT(eMesh)) THEN
         IF(ALLOCATED(eMesh)) DEALLOCATE(eMesh)
         ALLOCATE(eMesh(this%ne))

         DO ie = 1, this%ne
            eMesh(ie) = (ie-1)*del+eb
         ENDDO
      ENDIF

   END SUBROUTINE eMesh_gfinp

   SUBROUTINE eContour_gfinp(this,ef,irank,nz,e,de)

      USE m_constants
      USE m_grule
      USE m_ExpSave

      !Calculates the complex energy contour and
      !writes it into the corresponding arrays in gf

      CLASS(t_gfinp),       INTENT(IN)    :: this
      REAL,                 INTENT(IN)    :: ef
      INTEGER,              INTENT(IN)    :: irank
      INTEGER,              INTENT(IN)    :: nz
      COMPLEX, ALLOCATABLE, INTENT(INOUT) :: e(:)
      COMPLEX, ALLOCATABLE, INTENT(INOUT) :: de(:)

      INTEGER iz,n
      REAL e1, e2, sigma
      COMPLEX del
      REAL r, xr
      REAL, ALLOCATABLE :: x(:), w(:)

      !Make sure that the arrays are allocated
      IF(.NOT.(ALLOCATED(e).AND.ALLOCATED(de))) CALL juDFT_error("e or de array not allocated",calledby="eContour_gfinp")

      !Help arrays
      ALLOCATE(x(nz),source=0.0)
      ALLOCATE(w(nz),source=0.0)

      IF(this%mode.EQ.1) THEN

         sigma = this%sigma * pi_const
         IF(this%nmatsub > 0) THEN
            e1 = ef+this%eb
            n = 0

            !Left Vertical part (e1,0) -> (e1,sigma)
            del = this%nmatsub * CMPLX(0.0,sigma)
            CALL grule(this%n1,x(1:(this%n1)/2),w(1:(this%n1)/2))
            x = -x
            DO iz = 1, (this%n1+3)/2-1
               x(this%n1-iz+1) = -x(iz)
               w(this%n1-iz+1) =  w(iz)
            ENDDO
            DO iz = 1, this%n1
               n = n + 1
               IF(n.GT.nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               e(n)  = e1 + del + del * x(iz)
               de(n) = w(iz)*del
            ENDDO

            !Horizontal Part (eb,sigma) -> (et,sigma)
            del = (ef-30*this%sigma-e1)/2.0
            CALL grule(this%n2,x(1:(this%n2)/2),w(1:(this%n2)/2))
            x = -x
            DO iz = 1, (this%n2+3)/2-1
               x(this%n2-iz+1) = -x(iz)
               w(this%n2-iz+1) =  w(iz)
            ENDDO
            DO iz = 1, this%n2
               n = n + 1
               IF(n.GT.nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               e(n)  = del*x(iz) + del + e1 + 2 * this%nmatsub * ImagUnit * sigma
               de(n) = del*w(iz)
            ENDDO

            !Right Vertical part (et,sigma) -> infty
            CALL grule(this%n3,x(1:(this%n3)/2),w(1:(this%n3)/2))
            x = -x
            DO iz = 1, (this%n3+3)/2-1
               x(this%n3-iz+1) = -x(iz)
               w(this%n3-iz+1) =  w(iz)
            ENDDO
            del = 30*this%sigma
            DO iz = 1, this%n3
               n = n + 1
               IF(n.GT.nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               e(n)  = del*x(iz)+ef +  2 * this%nmatsub * ImagUnit * sigma
               de(n) = w(iz)*del/(1.0+exp_save((REAL(e(n))-ef)/this%sigma))
            ENDDO

            !Matsubara frequencies
            DO iz = this%nmatsub , 1, -1
               n = n + 1
               IF(n.GT.nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               e(n)  = ef + (2*iz-1) * ImagUnit *sigma
               de(n) =  -2 * ImagUnit * sigma
            ENDDO
         ENDIF
      ELSE IF(this%mode.EQ.2) THEN

         !Semicircle
         e1 = ef+this%eb
         e2 = ef+this%et

         !Radius
         r  = (e2-e1)*0.5
         !midpoint
         xr = (e2+e1)*0.5

         CALL grule(this%ncirc,x(1:(this%ncirc)/2),w(1:(this%ncirc)/2))

         DO iz = 1, this%ncirc/2
            x(this%ncirc-iz+1) = -x(iz)
            w(this%ncirc-iz+1) =  w(iz)
         ENDDO
         DO iz = 1, this%ncirc
            e(iz)  = xr + ImagUnit * r * EXP(ImagUnit*pi_const/2.0 * x(iz))
            de(iz) = pi_const/2.0 * r * w(iz) * EXP(ImagUnit*pi_const/2.0 * x(iz))
            !Scale the imaginary part with the given factor alpha
            e(iz)  = REAL(e(iz))  + ImagUnit * this%alpha * AIMAG(e(iz))
            de(iz) = REAL(de(iz)) + ImagUnit * this%alpha * AIMAG(de(iz))
         ENDDO

      ELSE IF(this%mode.EQ.3) THEN

         !Equidistant contour (without vertical edges)
         del = (this%et-this%eb)/REAL(nz-1)
         DO iz = 1, this%nDOS
            e(iz) = (iz-1) * del + this%eb + ImagUnit * this%sigmaDOS
            IF(this%l_dosfermi) THEN
               de(iz) = del * 1.0/(1.0+exp_save((REAL(e(iz))-ef)/this%sigmaDOS))
            ELSE
               de(iz) = del
            ENDIF
         ENDDO

      ELSE
         CALL juDFT_error("Invalid mode for energy contour in Green's function calculation", calledby="eContour_gfinp")
      END IF

      IF(irank.EQ.0) THEN
         !Write out the information about the energy contour
         WRITE(6,"(A)") "---------------------------------------------"
         WRITE(6,"(A)") " Green's function energy contour"
         WRITE(6,"(A)") "---------------------------------------------"
         WRITE(6,1000) this%mode
         WRITE(6,*)

         SELECT CASE(this%mode)

         CASE(1)
            WRITE(6,"(A)") "Rectangular Contour: "
            WRITE(6,1010) nz, this%nmatsub,this%n1,this%n2,this%n3
            WRITE(6,"(A)") "Energy limits (rel. to fermi energy): "
            WRITE(6,1040) this%eb,0.0
         CASE(2)
            WRITE(6,"(A)") "Semicircle Contour: "
            WRITE(6,1020) nz, this%alpha
            WRITE(6,"(A)") "Energy limits (rel. to fermi energy): "
            WRITE(6,1040) this%eb,this%et
         CASE(3)
            WRITE(6,"(A)") "Equidistant Contour for DOS calculations: "
            WRITE(6,1030) nz, this%sigmaDOS
            WRITE(6,"(A)") "Energy limits (rel. to fermi energy): "
            WRITE(6,1040) this%eb,this%et
         CASE default

         END SELECT

         !Write out points and weights
         WRITE(6,*)
         WRITE(6,"(A)") " Energy points: "
         WRITE(6,"(A)") "---------------------------------------------"
         DO iz = 1, nz
            WRITE(6,1050) REAL(e(iz)), AIMAG(e(iz)), REAL(de(iz)), AIMAG(de(iz))
         ENDDO

1000     FORMAT("Using energy contour mode: ", I1)
1010     FORMAT("nz: ", I5.1,"; nmatsub: ", I5.1,"; n1: ", I5.1,"; n2: ", I5.1,"; n3: ", I5.1)
1020     FORMAT("nz: ", I5.1," alpha: ", f8.4)
1030     FORMAT("n: ", I5.1,"; sigma: ", f8.4)
1040     FORMAT("eb: ", f8.4,"; et: ",f8.4)
1050     FORMAT(2f8.4,"      weight: ",2e15.4)
      ENDIF

   END SUBROUTINE eContour_gfinp

END MODULE m_types_gfinp
