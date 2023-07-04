!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts
   USE m_judft
   USE m_types_fleurinput_base
   USE m_constants
   IMPLICIT NONE
   PRIVATE
   type t_eibz
      integer :: nkpt
      integer, allocatable :: pointer(:)
   contains
      PROCEDURE :: init => init_eibz
      procedure :: calc_nkpt_EIBZ    => calc_nkpt_EIBZ
      procedure :: calc_pointer_EIBZ => calc_pointer_EIBZ
   end type t_eibz

   TYPE, EXTENDS(t_fleurinput_base):: t_kpts
      INTEGER                        :: kptsKind = 0
      CHARACTER(len=40)              :: kptsName = "default"
      character(len=100)             :: comment = ""
      INTEGER                        :: nkpt = 0
      INTEGER                        :: nkpt3(3) = 0 ! This variable is not supposed to be reliable. It is only meant as additional user information with the available data at IO time.
      INTEGER                        :: ntet = 0
      LOGICAL                        :: l_gamma = .FALSE.
      !(3,nkpt) k-vectors internal units
      REAL, ALLOCATABLE              :: bk(:, :) !(xyz,nk)
      !(nkpts) weights
      REAL, ALLOCATABLE              :: wtkpt(:)
      INTEGER                        :: nkptf = 0   !<k-vectors in full BZ
      REAL, ALLOCATABLE              :: bkf(:, :) !(xyz,nk)
      INTEGER, ALLOCATABLE           :: bkp(:)
      INTEGER, ALLOCATABLE           :: bksym(:)
      INTEGER                        :: numSpecialPoints = 0
      INTEGER, ALLOCATABLE           :: specialPointIndices(:)
      CHARACTER(LEN=50), ALLOCATABLE :: specialPointNames(:)
      REAL, ALLOCATABLE              :: specialPoints(:, :)
      INTEGER, ALLOCATABLE           :: ntetra(:, :)
      INTEGER, ALLOCATABLE           :: tetraList(:,:) !List with the tetrahedra containing the current kpoint (for more efficient loops)
      REAL, ALLOCATABLE              :: voltet(:)
      REAL, ALLOCATABLE              :: sc_list(:, :) !list for all information about folding of bandstructure (need for unfoldBandKPTS)((k(x,y,z),K(x,y,z),m(g1,g2,g3)),(nkpt),k_original(x,y,z))
      type(t_eibz), allocatable      :: EIBZ(:)
      logical                        :: l_set_eibz = .False.
   CONTAINS
      PROCEDURE :: calcCommonFractions
      PROCEDURE :: add_special_line
      PROCEDURE :: print_xml
      PROCEDURE :: read_xml_kptsByIndex
      PROCEDURE :: read_kpts_by_name
      PROCEDURE :: read_xml => read_xml_kpts
      PROCEDURE :: mpi_bc => mpi_bc_kpts
      procedure :: get_nk => kpts_get_nk
      procedure :: to_first_bz => kpts_to_first_bz
      procedure :: is_kpt => kpts_is_kpt
      procedure :: init => init_kpts
      procedure :: initTetra
      PROCEDURE :: tetrahedron_regular
      procedure :: calcNkpt3 => nkpt3_kpts
      procedure :: find_gamma => find_gamma_kpts
   ENDTYPE t_kpts

   PUBLIC :: t_kpts
CONTAINS

   function kpts_get_nk(kpts, kpoint) result(ret_idx)
      ! get the index of a kpoint
      implicit NONE
      class(t_kpts), intent(in)    :: kpts
      real, intent(in)             :: kpoint(3)
      integer                      :: idx, ret_idx
      real                         :: kpt_bz(3)

      kpt_bz = kpts%to_first_bz(kpoint)

      ret_idx = 0
      DO idx = 1, kpts%nkptf
         IF (all(abs(kpt_bz - kpts%bkf(:, idx)) < 1E-06)) THEN
            ret_idx = idx
            exit
         END IF
      END DO
   end function kpts_get_nk

   function kpts_to_first_bz(kpts, k_in) result(k)
      implicit NONE
      class(t_kpts), intent(in)  :: kpts
      real, intent(in)           :: k_in(3)
      real                       :: k(3)

      k = k_in
      ! everything close to 0 or 1 get's mapped to 0 and 1
      where (abs(k - dnint(k)) < 1e-8) k = dnint(k)

      ! map to 0 -> 1 interval
      k = k - floor(k)
   end function kpts_to_first_bz

   function kpts_is_kpt(kpts, kpoint) result(is_kpt)
      implicit none
      class(t_kpts), intent(in)  :: kpts
      real, intent(in)           :: kpoint(3)
      logical                    :: is_kpt

      is_kpt = kpts%get_nk(kpoint) > 0
   end function kpts_is_kpt
   SUBROUTINE mpi_bc_kpts(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_kpts), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF

      CALL mpi_bc(this%nkpt, rank, mpi_comm)
      CALL mpi_bc(this%kptsKind, rank, mpi_comm)
      CALL mpi_bc(this%nkpt3(1), rank, mpi_comm)
      CALL mpi_bc(this%nkpt3(2), rank, mpi_comm)
      CALL mpi_bc(this%nkpt3(3), rank, mpi_comm)
      CALL mpi_bc(this%ntet, rank, mpi_comm)
      CALL mpi_bc(this%l_gamma, rank, mpi_comm)
      CALL mpi_bc(this%bk, rank, mpi_comm)
      CALL mpi_bc(this%wtkpt, rank, mpi_comm)
      CALL mpi_bc(this%nkptf, rank, mpi_comm)
      CALL mpi_bc(this%bkf, rank, mpi_comm)
      CALL mpi_bc(this%bkp, rank, mpi_comm)
      CALL mpi_bc(this%bksym, rank, mpi_comm)
      CALL mpi_bc(this%numSpecialPoints, rank, mpi_comm)
      CALL mpi_bc(this%specialPointIndices, rank, mpi_comm)
      CALL mpi_bc(this%specialPoints, rank, mpi_comm)
      CALL mpi_bc(this%ntetra, rank, mpi_comm)
      CALL mpi_bc(this%tetraList,rank,mpi_comm)
      CALL mpi_bc(this%voltet, rank, mpi_comm)
      CALL mpi_bc(this%sc_list, rank, mpi_comm)
   END SUBROUTINE mpi_bc_kpts

   recursive logical function read_kpts_by_name(this,filename,name)
      USE m_calculator
      CLASS(t_kpts), INTENT(inout):: this
      character(len=*),INTENT(IN) :: filename,name
      
      character(len=150):: line
      integer           :: error,n,fid

      OPEN(newunit=fid,file=filename,action='READ')
      read_kpts_by_name=.false.
      DO while(.not. read_kpts_by_name)
         read(fid,"(a)",iostat=error) line
         IF (error.ne.0) exit 
         IF (index(line,"kPointList ")>0) THEN
            line=line(index(line,'name="')+6:)
            line=line(:index(line,'"')-1)
            if (line==trim(name)) THEN
               !Found kpointlist with correct name
               DO n=1,this%nkpt
                  read(fid,"(a)") line
                  line=line(index(line,"weight"):)
                  line=line(index(line,'"')+1:)
                  this%wtkpt(n)=evaluateFirstOnly(line(:index(line,'"')-1))
                  line=line(index(line,'>')+1:index(line,'<')-1)
                  this%bk(1, n) = evaluatefirst(line)
                  this%bk(2, n) = evaluatefirst(line)
                  this%bk(3, n) = evaluatefirst(line)
               ENDDO   
               read_kpts_by_name=.true.
            endif   
         ENDIF
         if (index(line,'<xi:include xmlns:xi="http://www.w3.org/2001/XInclude"')>0) THEN
            line=line(index(line,'href="')+6:)
            line=line(:index(line,'"')-1)
            read_kpts_by_name=this%read_kpts_by_name(line,name)
         endif
      enddo
      close(fid)
   end function      

   SUBROUTINE read_xml_kptsByIndex(this, xml, kptsIndex)
      USE m_types_xml
      USE m_calculator
      CLASS(t_kpts), INTENT(inout):: this
      TYPE(t_xml), INTENT(INOUT) ::xml
      INTEGER, INTENT(IN), OPTIONAL :: kptsIndex

      INTEGER :: numNodes, i, number_sets, n, currentIndex
      LOGICAL :: foundList, l_band
      CHARACTER(len=200)::str, path, path2, label, altPurpose
      CHARACTER(LEN=40) :: listName, typeString
      CHARACTER(LEN=255),ALLOCATABLE :: tetra_string(:)
      IF (xml%versionNumber > 31) then
        WRITE (path, "(a,i0,a)") '/fleurInput/cell/bzIntegration/kPointLists/kPointList[', kptsIndex, ']'
      ELSE
        path='/fleurInput/calculationSetup/bzIntegration/kPointList'
      endif
      numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path)))
      IF(numNodes.NE.1) THEN
         WRITE(*,*) 'kPointList index is ', kptsIndex
         CALL judft_error(("kPointList for index is not available."))
      END IF
      IF (xml%versionNumber > 31) THEN
         listName = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(path)//'/@name')))
         this%kptsName = TRIM(ADJUSTL(listName))

         this%kptsKind = 0
         this%nkpt3(:) = 0
         typeString = xml%GetAttributeValue(TRIM(path)//'/@type')
         SELECT CASE(TRIM(ADJUSTL(typeString)))
            CASE ('unspecified')
               this%kptsKind = KPTS_KIND_UNSPECIFIED
            CASE ('mesh')
               this%kptsKind = KPTS_KIND_MESH
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@nx')
               IF(numNodes.EQ.1) this%nkpt3(1) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@nx'))
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@ny')
               IF(numNodes.EQ.1) this%nkpt3(2) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@ny'))
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@nz')
               IF(numNodes.EQ.1) this%nkpt3(3) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@nz'))
            CASE ('path')
               this%kptsKind = KPTS_KIND_PATH
            CASE ('tria')
               this%kptsKind = KPTS_KIND_TRIA
            CASE ('tria-bulk')
               this%kptsKind = KPTS_KIND_TRIA_BULK
            CASE ('SPEX-mesh') ! (this is deprecated)
               this%kptsKind = KPTS_KIND_SPEX_MESH
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@nx')
               IF(numNodes.EQ.1) this%nkpt3(1) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@nx'))
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@ny')
               IF(numNodes.EQ.1) this%nkpt3(2) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@ny'))
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@nz')
               IF(numNodes.EQ.1) this%nkpt3(3) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@nz'))
            CASE ('spex-mesh') ! (same as 'SPEX-mesh')
               this%kptsKind = KPTS_KIND_SPEX_MESH
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@nx')
               IF(numNodes.EQ.1) this%nkpt3(1) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@nx'))
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@ny')
               IF(numNodes.EQ.1) this%nkpt3(2) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@ny'))
               numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/@nz')
               IF(numNodes.EQ.1) this%nkpt3(3) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(path)//'/@nz'))
            CASE DEFAULT
               this%kptsKind = KPTS_KIND_UNSPECIFIED
               WRITE(*,*) 'WARNING: Unknown k point list type. Assuming "unspecified"'
         END SELECT
      END IF
      this%nkpt = evaluateFirstOnly(xml%GetAttributeValue(TRIM(path)//'/@count'))
      numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/kPoint')
      IF (numNodes.NE.this%nkpt) THEN
         CALL judft_error("Inconsistent number of k-points in kPointList: "//TRIM(ADJUSTL(this%kptsName)))
      END IF

      ! count special points
      this%numSpecialPoints = 0
      IF (this%kptskind==KPTS_KIND_PATH) THEN
         DO i = 1, numNodes
            label = ''
            WRITE (path2, "(a,a,i0,a)") TRIM(ADJUSTL(path)), "/kPoint[", i, "]"
            label = xml%GetAttributeValue(TRIM(path2)//'/@label')
            IF (TRIM(ADJUSTL(label)).NE.'') this%numSpecialPoints = this%numSpecialPoints + 1
         END DO
      ENDIF   

      IF(this%numSpecialPoints.NE.0) THEN
         IF(ALLOCATED(this%specialPointIndices)) DEALLOCATE(this%specialPointIndices)
         ALLOCATE(this%specialPointIndices(this%numSpecialPoints))
         ALLOCATE (this%specialPoints(3, this%numSpecialPoints))
         ALLOCATE (this%specialPointNames(this%numSpecialPoints))
         ! set special points
         currentIndex = 0
         DO i = 1, numNodes
            label = ''
            WRITE (path2, "(a,a,i0,a)") TRIM(ADJUSTL(path)), "/kPoint[", i, "]"
            label = xml%GetAttributeValue(TRIM(path2)//'/@label')
            IF (TRIM(ADJUSTL(label)).NE.'') THEN
               currentIndex = currentIndex + 1
               this%specialPointNames(currentIndex) = TRIM(ADJUSTL(label))
               this%specialPointIndices(currentIndex) = i
               str = xml%getAttributeValue(TRIM(ADJUSTL(path2)))
               this%specialPoints(1, currentIndex) = evaluatefirst(str)
               this%specialPoints(2, currentIndex) = evaluatefirst(str)
               this%specialPoints(3, currentIndex) = evaluatefirst(str)
            END IF
         END DO
      END IF

      ALLOCATE (this%bk(3, this%nkpt))
      ALLOCATE (this%wtkpt(this%nkpt))
      if (.not. this%read_kpts_by_name("inp.xml",this%kptsName)) THEN
         print *,"WARNING, new k-point reader could not be used. Please check your inp.xml/kpts.xml"
         DO i = 1, this%nkpt
            WRITE (path2, "(a,a,i0,a)") TRIM(ADJUSTL(path)), "/kPoint[", i, "]"
            this%wtkpt(i) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(path2))//'/@weight',.true.))
            str = xml%getAttributeValue(TRIM(ADJUSTL(path2)),.true.)
            this%bk(1, i) = evaluatefirst(str)
            this%bk(2, i) = evaluatefirst(str)
            this%bk(3, i) = evaluatefirst(str)
         END DO
      endif   

!      n = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/tetraeder')
!      IF (n .EQ. 1) THEN
!         this%ntet = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/tetraeder/tet')
!         ALLOCATE(tetra_string(this%ntet))
!         call xml%GetAttributeValue_List(TRIM(ADJUSTL(path))//'/tetraeder/tet',tetra_string)
!         ALLOCATE (this%voltet(this%ntet), this%ntetra(4, this%ntet))
!         DO n = 1, this%ntet
!            WRITE (path2, "(a,a,i0,a)") TRIM(ADJUSTL(path)), "/tetraeder/tet[", n, "]"
!            this%voltet(n) = Evaluatefirst(Tetra_string(N))
!            this%ntetra(1,n) = Evaluatefirst(Tetra_string(N))
!            this%ntetra(2,n) = Evaluatefirst(Tetra_string(N))
!            this%ntetra(3,n) = Evaluatefirst(Tetra_string(N))
!            this%ntetra(4,n) = Evaluatefirst(Tetra_string(N))
!
!                        !str = xml%getAttributeValue(TRIM(ADJUSTL(path2)),.true.)
!            !READ (str,*) this%ntetra(:,n)
!         ENDDO
!         deallocate(tetra_string)
!      ENDIF
!
!      n = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/triangles')
!      IF (n .EQ. 1) THEN
!         this%ntet = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/triangles/tria')
!         ALLOCATE(tetra_string(this%ntet))
!         call xml%GetAttributeValue_List(TRIM(ADJUSTL(path))//'/triangles/tria',tetra_string)
!         ALLOCATE (this%voltet(this%ntet), this%ntetra(3, this%ntet))
!         DO n = 1, this%ntet
!            this%voltet(n) = evaluateFirst(Tetra_string(n))
!            this%ntetra(1,n) = Evaluatefirst(Tetra_string(N))
!            this%ntetra(2,n) = Evaluatefirst(Tetra_string(N))
!            this%ntetra(3,n) = Evaluatefirst(Tetra_string(N))
!         ENDDO
!      ENDIF
      this%wtkpt = this%wtkpt/sum(this%wtkpt) !Normalize k-point weight
   END SUBROUTINE read_xml_kptsByIndex

   SUBROUTINE read_xml_kpts(this, xml)
      USE m_types_xml
      USE m_calculator
      CLASS(t_kpts), INTENT(inout):: this
      TYPE(t_xml), INTENT(INOUT) ::xml

      INTEGER :: numNodes, i, number_sets, n, currentIndex, kptsIndex
      LOGICAL :: foundList, l_band
      CHARACTER(len=200)::str, path, path2, label, altPurpose
      CHARACTER(LEN=40) :: listName, typeString


      IF (xml%versionNumber > 31) then
        listName = xml%GetAttributeValue('/fleurInput/cell/bzIntegration/kPointListSelection/@listName')

        numNodes = xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/altKPointList')

        l_band = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@band'))
        IF (l_band) THEN
          DO i = 1 , numNodes
            WRITE (path2, "(a,i0,a)") '/fleurInput/cell/bzIntegration/altKPointList[',i,']'
            altPurpose = ''
            altPurpose = xml%GetAttributeValue(TRIM(ADJUSTL(path2))//'/@purpose')
            IF (TRIM(ADJUSTL(altPurpose)).EQ.'bands') THEN
              listName = ''
              listName = xml%GetAttributeValue(TRIM(ADJUSTL(path2))//'/@listName')
              EXIT
            END IF
          END DO
        END IF

        numNodes = xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/kPointLists/kPointList')
        foundList = .FALSE.
        kptsIndex = 0
        DO i = 1, numNodes
          path = ''
          WRITE (path, "(a,i0,a)") '/fleurInput/cell/bzIntegration/kPointLists/kPointList[', i, ']'
          IF (TRIM(ADJUSTL(listName)) == xml%GetAttributeValue(TRIM(path)//'/@name')) THEN
            kptsIndex = i
            foundList = .TRUE.
            EXIT
          END IF
        END DO

      ELSE
        foundList=.true.
        kptsIndex=1
      endif
      !simply read the first k-list

      IF(.NOT.foundList) THEN
        CALL judft_error(("No kPointList named: "//TRIM(ADJUSTL(listName))//" found."))
      END IF

      CALL read_xml_kptsByIndex(this, xml, kptsIndex)

      IF (l_band) THEN
         IF (.NOT.this%kptsKind.EQ.KPTS_KIND_PATH) THEN
            CALL juDFT_warn('Chosen k-point list is not compatible to band structure calculations.', calledby='read_xml_kpts')
         END IF
      ELSE
         IF (this%kptsKind.EQ.KPTS_KIND_PATH) THEN
            CALL juDFT_warn('Chosen k-point list is only compatible to band structure calculations.', calledby='read_xml_kpts')
         END IF
      END IF

   END SUBROUTINE read_xml_kpts

   SUBROUTINE print_xml(kpts, kptsUnit, filename)
      CLASS(t_kpts), INTENT(in):: kpts
      INTEGER, INTENT(in)         :: kptsUnit
      CHARACTER(len=*), INTENT(in), OPTIONAL::filename

      INTEGER :: n, iSpecialPoint, i, nkq_pairs
      REAL :: commonFractions(3)
      LOGICAL :: l_exist
      CHARACTER(LEN=20) :: posString(3)
      CHARACTER(LEN=50) :: label

      label = ''

      IF(.NOT.ALLOCATED(kpts%bk)) RETURN

      IF (PRESENT(filename)) THEN
         INQUIRE (file=filename, exist=l_exist)
         IF (l_exist) THEN
            OPEN (kptsUnit, file=filename, action="write", position="append")
         ELSE
            OPEN (kptsUnit, file=filename, action="write")
         END IF
      ENDIF

      commonFractions(:) = -1.0

205   FORMAT('            <kPointList name="', a, '" count="', i0, '" type="', a, '">')
2051  FORMAT('            <kPointList name="', a, '" count="', i0, '" nx="', i0, '" ny="', i0, '" nz="', i0,  '" type="', a, '">')
2052  FORMAT('            <kPointList name="', a, '" count="', i0, '" nx="', i0, '" ny="', i0, '" nz="', i0, &
                                                                                          '" nkq_pairs="', i0, '" type="', a, '">')
      IF(kpts%kptsKind.EQ.KPTS_KIND_MESH) THEN
         if(kpts%l_gamma .and. kpts%l_set_eibz) then 
            nkq_pairs = sum([(kpts%eibz(i)%nkpt, i=1, size(kpts%eibz))])
            WRITE (kptsUnit, 2052) TRIM(ADJUSTL(kpts%kptsName)), kpts%nkpt, kpts%nkpt3(1), kpts%nkpt3(2), kpts%nkpt3(3),&
                                                                        nkq_pairs, TRIM(ADJUSTL(kptsKindString_consts(kpts%kptsKind)))
         else
            WRITE (kptsUnit, 2051) TRIM(ADJUSTL(kpts%kptsName)), kpts%nkpt, kpts%nkpt3(1), kpts%nkpt3(2), kpts%nkpt3(3), TRIM(ADJUSTL(kptsKindString_consts(kpts%kptsKind)))
         endif
         CALL calcCommonFractions(kpts,commonFractions)
      ELSE
         WRITE (kptsUnit, 205) TRIM(ADJUSTL(kpts%kptsName)), kpts%nkpt, TRIM(ADJUSTL(kptsKindString_consts(kpts%kptsKind)))
      END IF

      DO n = 1, kpts%nkpt
         DO iSpecialPoint = 1, kpts%numSpecialPoints
            IF (kpts%specialPointIndices(iSpecialPoint).EQ.n) THEN
               label = kpts%specialPointNames(iSpecialPoint)
               EXIT
            END IF
         END DO

         posString(:) = ''
         IF((kpts%kptsKind.EQ.KPTS_KIND_MESH).AND.(ALL(commonFractions(:).GT.1.0))) THEN
            WRITE(posString(1),'(f7.2,a,f0.2)') commonFractions(1)*kpts%bk(1, n), '/' , commonFractions(1)
            WRITE(posString(2),'(f7.2,a,f0.2)') commonFractions(2)*kpts%bk(2, n), '/' , commonFractions(2)
            WRITE(posString(3),'(f7.2,a,f0.2)') commonFractions(3)*kpts%bk(3, n), '/' , commonFractions(3)
         ELSE
            WRITE(posString(1),'(f19.16)') kpts%bk(1, n)
            WRITE(posString(2),'(f19.16)') kpts%bk(2, n)
            WRITE(posString(3),'(f19.16)') kpts%bk(3, n)
         END IF

206      FORMAT('               <kPoint weight="', f20.13, '">', a, ' ', a, ' ', a, '</kPoint>')
2061     FORMAT('               <kPoint weight="', f20.13, '" label="',a ,'">', a, ' ', a, ' ', a, '</kPoint>')
         IF(label.EQ.'') THEN
            WRITE (kptsUnit, 206) kpts%wtkpt(n), TRIM(posString(1)), TRIM(posString(2)), TRIM(posString(3))
         ELSE
            WRITE (kptsUnit, 2061) kpts%wtkpt(n), TRIM(ADJUSTL(label)), TRIM(posString(1)), TRIM(posString(2)), TRIM(posString(3))
         END IF
         label = ''
      END DO
!      IF (kpts%ntet > 0) THEN
!         IF (SIZE(kpts%ntetra, 1).EQ.3) THEN
!            !Film --> Triangles
!            WRITE (kptsUnit, 209) kpts%ntet
!209         FORMAT('               <triangles ntria="', i0, '">')
!            DO n = 1, kpts%ntet
!210            FORMAT('                  <tria>', f20.13, i0, ' ', i0, ' ', i0, '</tria>')
!               WRITE (kptsUnit, 210) kpts%voltet(n), kpts%ntetra(:, n)
!            END DO
!            WRITE (kptsUnit, '(a)') '               </triangles>'
!         ENDIF
!      ELSE
!         DO n = 1, kpts%numSpecialPoints
!            WRITE (kptsUnit, 211) TRIM(ADJUSTL(kpts%specialPointNames(n))), kpts%specialPoints(:, n)
!211         FORMAT('            <specialPoint name="', a, '">', f10.6, ' ', f10.6, ' ', f10.6, '</specialPoint>')
!         END DO
!      END IF
      WRITE (kptsUnit, '(a)') ('            </kPointList>')
      IF (PRESENT(filename)) CLOSE (kptsUnit)
   END SUBROUTINE print_xml

   SUBROUTINE add_special_line(kpts, point, name)
      CLASS(t_kpts), INTENT(inout):: kpts
      CHARACTER(len=*), INTENT(in):: name
      REAL, INTENT(in)            :: point(3)

      CHARACTER(len=50), ALLOCATABLE:: names(:)
      REAL, ALLOCATABLE             :: points(:, :)

      IF (kpts%numspecialpoints > 0) THEN
         CALL MOVE_ALLOC(kpts%specialPointNames, names)
         CALL MOVE_ALLOC(kpts%specialPoints, points)
         DEALLOCATE (kpts%specialpointindices)
         ALLOCATE (kpts%specialPoints(3, SIZE(names) + 1))
         ALLOCATE (kpts%specialPointIndices(SIZE(names) + 1))
         ALLOCATE (kpts%specialPointNames(SIZE(names) + 1))
         kpts%specialPoints(:, :SIZE(names)) = points
         kpts%specialPointNames(:SIZE(names)) = names
      ELSE
         ALLOCATE (kpts%specialPoints(3, 1))
         ALLOCATE (kpts%specialPointIndices(1))
         ALLOCATE (kpts%specialPointNames(1))
      ENDIF
      kpts%numspecialpoints = kpts%numspecialpoints + 1
      kpts%specialPoints(:, kpts%numspecialpoints) = point
      kpts%specialPointNames(kpts%numspecialpoints) = name
   END SUBROUTINE add_special_line

   subroutine init_EIBZ(EIBZ, kpts, sym, nk)
      USE m_types_sym
      implicit none
      class(t_eibz), intent(inout)   :: EIBZ
      type(t_kpts), intent(in)       :: kpts
      type(t_sym), intent(in)        :: sym
      integer, intent(in)            :: nk

      call eibz%calc_nkpt_EIBZ(kpts, sym, nk)
      call eibz%calc_pointer_EIBZ(kpts, sym, nk)
   end subroutine init_EIBZ

   SUBROUTINE initTetra(kpts,input,cell,sym,l_soc_or_ss)
      USE m_juDFT
      USE m_constants
      USE m_types_input
      USE m_types_cell
      USE m_types_sym
      USE m_tetcon
      USE m_triang
      CLASS(t_kpts),    INTENT(INOUT) :: kpts
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_sym),      INTENT(IN)    :: sym
      LOGICAL,          INTENT(IN)    :: l_soc_or_ss

      INTEGER, PARAMETER :: nop48  = 48

      INTEGER :: i, j, ikpt, ntet, itet
      INTEGER :: ndiv3,nsym,addSym
      REAL    :: volirbz, as
      REAL    :: vkxyz(3,kpts%nkpt)
      LOGICAL :: l_tria

      REAL    :: bltv(3,3)          ! cartesian Bravais lattice basis (a.u.)
      REAL    :: rltv(3,3)          ! reciprocal lattice basis (2\pi/a.u.)
      REAL    :: ccr(3,3,nop48)     ! rotation matrices in cartesian repr.
      REAL    :: rlsymr(3,3,nop48)  ! rotation matrices in reciprocal lattice basis representation
      REAL    :: binv(3,3)

      INTEGER, ALLOCATABLE :: ntetra(:,:) ! corners of the tetrahedrons
      REAL,    ALLOCATABLE :: voltet(:)   ! voulmes of the tetrahedrons

      nsym = MERGE(sym%nop2,sym%nop,input%film)
      bltv=TRANSPOSE(cell%amat)
      binv=TRANSPOSE(cell%bmat)/tpi_const

      DO i = 1, nsym
         rlsymr(:,:,i)=REAL(sym%mrot(:,:,i))
         ccr(:,:,i) = TRANSPOSE(MATMUL(MATMUL(binv(:,:),TRANSPOSE(rlsymr(:,:,i))),bltv(:,:)))
      END DO

      IF ((.NOT.l_soc_or_ss).AND.(2*nsym<nop48)) THEN
         IF ((input%film.AND.(.NOT.sym%invs2)).OR.((.NOT.input%film).AND.(.NOT.sym%invs))) THEN
            addSym = 0
            ! Note: We have to add the negative of each symmetry operation
            !       to exploit time reversal symmetry. However, if the new
            !       symmetry operation is the identity matrix it is excluded.
            !       This is the case iff it is (-Id) + a translation vector.
            DO i = 1, nsym
               ! This test assumes that ccr(:,:,1) is the identity matrix.
               IF(.NOT.ALL(ABS(ccr(:,:,1)+ccr(:,:,i)).LT.10e-10) ) THEN
                  ccr(:,:,nsym+addSym+1 ) = -ccr(:,:,i)
                  rlsymr(:,:,nsym+addSym+1 ) = -rlsymr(:,:,i)
                  addSym = addSym + 1
               END IF
            END DO
            nsym = nsym + addSym
         END IF
      END IF

      IF ((input%bz_integration.EQ.BZINT_METHOD_TRIA).AND.(.NOT.input%film)) THEN

         IF(kpts%kptsKind.NE.KPTS_KIND_TRIA_BULK) THEN
            CALL juDFT_error("'tria' tetrahedron decomposition for bulk systems needs a tria-bulk k-point set",&
                             calledby="initTetra")
         END IF

         DO j=1,kpts%nkpt
            vkxyz(:,j)=MATMUL(kpts%bk(:,j),cell%bmat)
         END DO
         ndiv3 = 6*(kpts%nkpt+1)

         ALLOCATE (ntetra(4,ndiv3))
         ALLOCATE (voltet(ndiv3))

         CALL tetcon(kpts%nkpt,ndiv3,cell%omtil,vkxyz,nsym, kpts%ntet,voltet,ntetra)

         WRITE (oUnit,'('' the number of tetrahedra '')')
         WRITE (oUnit,*) kpts%ntet
         WRITE (oUnit,'('' volumes of the tetrahedra '')')
         WRITE (oUnit,'(e19.12,1x,i5,5x,''voltet(i),i'')') (voltet(i),i,i=1,kpts%ntet)
         WRITE (oUnit,'('' corners of the tetrahedra '')')
         WRITE (oUnit, '(4(3x,4i4))') ((ntetra(j,i),j=1,4),i=1,kpts%ntet)
         WRITE (oUnit,'('' the # of different k-points '')')
         WRITE (oUnit,*) kpts%nkpt
         WRITE (oUnit,'('' k-points used to construct tetrahedra'')')
         WRITE (oUnit,'(3(4x,f10.6))') ((vkxyz(i,j),i=1,3),j=1,kpts%nkpt)

         volirbz =  tpi_const**3 /(real(nsym)*cell%omtil)
         DO i = 1, kpts%ntet
            voltet(i) = kpts%ntet * voltet(i) / volirbz 
         END DO

         IF(ALLOCATED(kpts%ntetra)) DEALLOCATE(kpts%ntetra)
         IF(ALLOCATED(kpts%voltet)) DEALLOCATE(kpts%voltet)
         ALLOCATE(kpts%ntetra(4,kpts%ntet))
         ALLOCATE(kpts%voltet(kpts%ntet))
         DO j = 1, kpts%ntet
            kpts%ntetra(1:4,j) = ntetra(1:4,j)
            kpts%voltet(j) = ABS(voltet(j))
         END DO
      END IF

      IF(input%bz_integration==BZINT_METHOD_TRIA .AND. input%film) THEN

         IF(kpts%kptsKind.NE.KPTS_KIND_MESH) THEN
            CALL juDFT_error("'tria' tetrahedron decomposition for film systems needs a k-point mesh",&
                             calledby="initTetra")
         END IF

         ALLOCATE (voltet(2*kpts%nkpt),ntetra(3,2*kpts%nkpt))
         voltet=0.0
         l_tria = .FALSE.
         CALL triang(kpts%bk,kpts%nkpt,ntetra,kpts%ntet,voltet,as,l_tria)
         !IF (sym%invs) THEN
         !   IF (abs(sym%nop2*as-0.5).GT.0.000001) l_tria=.false.
         !ELSE
         !   IF (abs(sym%nop2*as-1.0).GT.0.000001) l_tria=.false.
         !ENDIF
         !write(*,*) as,sym%nop2,l_tria

         !Match normalisation of other methods
         voltet = voltet/as*kpts%ntet

         IF(ALLOCATED(kpts%ntetra)) DEALLOCATE(kpts%ntetra)
         IF(ALLOCATED(kpts%voltet)) DEALLOCATE(kpts%voltet)
         ALLOCATE(kpts%ntetra(3,kpts%ntet))
         ALLOCATE(kpts%voltet(kpts%ntet))
         DO j = 1, kpts%ntet
            kpts%ntetra(1:3,j) = ntetra(1:3,j)
            kpts%voltet(j) = ABS(voltet(j))
         END DO
      END IF

      IF(input%bz_integration.EQ.BZINT_METHOD_TETRA) THEN
         !Regular decomposition of the Monkhorst Pack Grid into tetrahedra
         ! (kpts%init is supposed to be called before this point)
         IF((kpts%kptsKind.NE.KPTS_KIND_MESH).OR.(.NOT.kpts%l_gamma)) THEN
            CALL juDFT_error("Regular tetrahedron decomposition needs a gamma centered kpoint grid",&
                             calledby="initTetra")
         END IF
         CALL timestart('Tetrahedron decomposition')
         CALL kpts%tetrahedron_regular(input%film,cell,kpts%nkpt3,ntetra,voltet)
         CALL timestop('Tetrahedron decomposition')

         IF(ALLOCATED(kpts%ntetra)) DEALLOCATE(kpts%ntetra)
         IF(ALLOCATED(kpts%voltet)) DEALLOCATE(kpts%voltet)
         IF (.NOT.input%film) THEN
            ALLOCATE(kpts%ntetra(4,kpts%ntet))
            ALLOCATE(kpts%voltet(kpts%ntet))
            DO j = 1, kpts%ntet
               kpts%ntetra(:,j) = ntetra(1:4,j)
               kpts%voltet(j) = ABS(voltet(j))
            END DO
         ELSE
            ALLOCATE(kpts%ntetra(3,kpts%ntet))
            ALLOCATE(kpts%voltet(kpts%ntet))
            DO j = 1, kpts%ntet
               kpts%ntetra(:,j) = ntetra(1:3,j)
               kpts%voltet(j) = ABS(voltet(j))
            END DO
         END IF
      END IF

      IF((input%bz_integration.EQ.BZINT_METHOD_TETRA).OR.&
         (input%bz_integration.EQ.BZINT_METHOD_TRIA)) THEN
         CALL timestart("setup tetraList")
         allocate(kpts%tetraList( MERGE(2*sym%nop,sym%nop,.NOT.sym%invs)&
                                 *MERGE(6,24,input%film),kpts%nkpt),source=0)
         !$OMP parallel do default(none) private(ikpt,ntet,itet) shared(kpts)
         do ikpt = 1, kpts%nkpt
            ntet = 0
            do itet = 1, kpts%ntet
               IF(ANY(kpts%ntetra(:,itet).EQ.ikpt))THEN
                  ntet = ntet + 1
                  kpts%tetraList(ntet,ikpt) = itet
               ENDIF
            enddo
         enddo
         !$OMP end parallel do
         CALL timestop("setup tetraList")
      END IF

   END SUBROUTINE initTetra

   SUBROUTINE tetrahedron_regular(kpts,film,cell,grid,ntetra,voltet)
      USE m_types_cell
      USE m_juDFT
      USE m_constants
      CLASS(t_kpts),          INTENT(INOUT)  :: kpts
      LOGICAL,                INTENT(IN)     :: film
      TYPE(t_cell),           INTENT(IN)     :: cell
      INTEGER,                INTENT(IN)     :: grid(:)
      INTEGER, ALLOCATABLE,   INTENT(INOUT)  :: ntetra(:,:)
      REAL,    ALLOCATABLE,   INTENT(INOUT)  :: voltet(:)


      INTEGER :: ntetraCube,k1,k2,k3,ikpt,itetra,i,j,k,l
      INTEGER :: jtet,icorn,startIndex,itet
      REAL    :: vol,volbz,diag(2),minKpt(3)
      INTEGER :: iarr(3)
      LOGICAL :: l_new,l_equal_kpoints
      INTEGER, ALLOCATABLE :: tetra(:,:)
      INTEGER, ALLOCATABLE :: ntetraAll(:,:)
      INTEGER, ALLOCATABLE :: kcorn(:), kpt_tetra(:), ind(:)
      INTEGER, ALLOCATABLE :: p(:,:,:)

      !Determine the decomposition of each individual cube
      ! and the total volume of the brillouin zone
      IF(film) THEN
         volbz = cell%bmat(1,1)*cell%bmat(2,2)-cell%bmat(1,2)*cell%bmat(2,1)
         ALLOCATE(tetra(3,2),source=0)
         ALLOCATE(kcorn(4),source=0)
         ALLOCATE(kpt_tetra(3),source=0)
         ALLOCATE(ind(3),source=0)
         ntetraCube = 2
         tetra = reshape ( [ 1,2,3, 2,3,4], [ 3,2 ] )
         diag = cell%bmat(:2,2)/grid(:2) - cell%bmat(:2,1) / grid(:2)
         vol =  sum(diag*diag)/4.0
      ELSE
         volbz = ABS(det(cell%bmat))
         ALLOCATE(tetra(4,24),source=0)
         ALLOCATE(kcorn(8),source=0)
         ALLOCATE(kpt_tetra(4),source=0)
         ALLOCATE(ind(4),source=0)
         !Choose the tetrahedra decomposition along the shortest diagonal
         CALL get_tetra(cell,grid,ntetraCube,vol,tetra)
      ENDIF

      !We shift the k-points by this vector only for the pointer array
      DO i = 1, 3
         minKpt(i) = MINVAL(kpts%bkf(i,:))
      ENDDO

      !Set up pointer array for the kpts
      ALLOCATE(p(0:grid(1),0:grid(2),0:grid(3)),source=0)
      p = 0
      DO ikpt = 1, kpts%nkptf
         iarr = nint((kpts%bkf(:,ikpt)-minKpt(:))*grid)
         p(iarr(1),iarr(2),iarr(3)) = ikpt
      ENDDO
      p(grid(1),:,:) = p(0,:,:)
      p(:,grid(2),:) = p(:,0,:)
      p(:,:,grid(3)) = p(:,:,0)


      !Check for invalid indices
      IF(ANY(p<=0).OR.ANY(p>kpts%nkptf)) THEN
         CALL juDFT_error("Invalid kpoint index in pointer array",calledby="tetrahedron_regular")
      ENDIF

      !Temporary Size
      IF(film) THEN
         ALLOCATE(ntetra(3,kpts%nkptf*2),source=0)
         ALLOCATE(ntetraAll(3,kpts%nkptf*2),source=0)
         ALLOCATE(voltet(kpts%nkptf*2),source=0.0)
      ELSE
         ALLOCATE(ntetra(4,kpts%nkptf*6),source=0)
         ALLOCATE(ntetraAll(4,kpts%nkptf*6),source=0)
         ALLOCATE(voltet(kpts%nkptf*6),source=0.0)
      ENDIF

      !Set up the tetrahedrons
      !$omp parallel do default(none) &
      !$omp shared(grid,p,ntetraCube,kpts,tetra,film,ntetraAll) &
      !$omp private(k1,k2,k3,kcorn,itetra,startIndex,kpt_tetra,ind) &
      !$omp collapse(3)
      DO k3 = 0, MERGE(grid(3)-1,0,grid(3).NE.0)
         DO k2 = 0, grid(2)-1
            DO k1 = 0, grid(1)-1

               !Corners of the current cube
               kcorn(1) = p(k1  ,k2  ,k3  )
               kcorn(2) = p(k1+1,k2  ,k3  )
               kcorn(3) = p(k1  ,k2+1,k3  )
               kcorn(4) = p(k1+1,k2+1,k3  )
               IF(.NOT.film) THEN
                  kcorn(5) = p(k1  ,k2  ,k3+1)
                  kcorn(6) = p(k1+1,k2  ,k3+1)
                  kcorn(7) = p(k1  ,k2+1,k3+1)
                  kcorn(8) = p(k1+1,k2+1,k3+1)
               ENDIF

               !Now divide the cube into tetrahedra
               startIndex = (k3*grid(2)*grid(1)+k2*grid(1)+k1) * ntetraCube
               DO itetra = 1, ntetraCube
                  kpt_tetra = kpts%bkp(kcorn(tetra(:,itetra)))
                  ind = sort_int(kpts%bkp(kcorn(tetra(:,itetra))))
                  ntetraAll(:,startIndex + itetra) = kpt_tetra(ind(:))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$omp end parallel do

      !Check for symmetry equivalent tetrahedra
      kpts%ntet = 0
      DO itet = 1, ntetraCube*PRODUCT(grid(:MERGE(2,3,film)))
         l_new = .TRUE.
         tetraLoop: DO jtet = 1, kpts%ntet
            if(all(ntetraAll(:,itet)-ntetra(:,jtet).EQ.0)) then
               voltet(jtet) = voltet(jtet) + vol
               l_new = .false.
               exit tetraLoop
            endif
         ENDDO tetraLoop
         IF(l_new) THEN !This tetrahedron has no symmetry equivalents yet
            kpts%ntet = kpts%ntet+1
            ntetra(:,kpts%ntet) = ntetraAll(:,itet)
            voltet(kpts%ntet) = vol
         ENDIF
      ENDDO

      !Has the whole brillouin zone been covered?
      IF(ABS(SUM(voltet)-volbz).GT.1E-10) THEN
         CALL juDFT_error("tetrahedron_regular failed", calledby="tetrahedron_regular")
      ENDIF

      !Normalize volumes
      voltet = voltet/volbz

      !Rescale volumes for IO to inp.xml
      !(so weights dont get to small for IO with dense meshes)
      voltet = voltet * kpts%ntet
   
      contains

      pure function sort_int(arr) result(ind)

         integer, intent(in) :: arr(:)
         integer :: ind(size(arr))

         integer i,j
         integer tmp


         DO i = 1, size(arr)
            ind(i) = i
         ENDDO

         DO i = 1, size(arr)-1
            DO j = i+1, size(arr)
               IF (arr(ind(i)).GT.arr(ind(j))) THEN
                  tmp = ind(i)
                  ind(i) = ind(j)
                  ind(j) = tmp
               ENDIF
            ENDDO
         ENDDO

      end function

   END SUBROUTINE tetrahedron_regular

   SUBROUTINE get_tetra(cell,grid,ntetra,vol,tetra)
      USE m_types_cell
      USE m_juDFT
      USE m_constants
      TYPE(t_cell),  INTENT(IN)     :: cell
      INTEGER,       INTENT(IN)     :: grid(:)
      INTEGER,       INTENT(INOUT)  :: ntetra
      REAL,          INTENT(INOUT)  :: vol
      INTEGER,       INTENT(INOUT)  :: tetra(:,:)

      REAL rlv(3,3),diag(4),d(3)
      INTEGER idmin

      !Calculate the lengths of the three diagonals
      rlv(:,1) = cell%bmat(:,1) / grid
      rlv(:,2) = cell%bmat(:,2) / grid
      rlv(:,3) = cell%bmat(:,3) / grid

      vol = 1/6.0*ABS(det(rlv))
      d = rlv(:,1) + rlv(:,3) - rlv(:,2)
      diag(1) = sum(d*d)
      d = rlv(:,2) + rlv(:,3) - rlv(:,1)
      diag(2) = sum(d*d)
      d = rlv(:,1) + rlv(:,2) + rlv(:,3)
      diag(3) = sum(d*d)
      d = rlv(:,1) + rlv(:,2) - rlv(:,3)
      diag(4) = sum(d*d)
      idmin = minloc(diag,1)

      ntetra = 0
      !From spex tetrahedron.f (For now we only choose one decomposition)
      if(idmin==1) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 1,2,3,6, 5,7,3,6, 1,5,3,6, 2,4,3,6, 4,8,3,6, 7,8,3,6 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(idmin==2) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 5,6,2,7, 1,5,2,7, 1,3,2,7, 8,6,2,7, 4,3,2,7, 8,4,2,7 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(idmin==3) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 2,6,1,8, 2,4,1,8, 3,4,1,8, 3,7,1,8, 5,7,1,8, 5,6,1,8 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(idmin==4) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 2,6,4,5, 1,2,4,5, 1,3,4,5, 3,7,4,5, 7,8,4,5, 6,8,4,5 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif

   END SUBROUTINE get_tetra

   REAL FUNCTION det(m)
      REAL m(:,:)
      det = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
            m(2,1)*m(3,2)*m(1,3) - m(1,3)*m(2,2)*m(3,1) - &
            m(2,3)*m(3,2)*m(1,1) - m(2,1)*m(1,2)*m(3,3)
   END FUNCTION det


   SUBROUTINE init_kpts(kpts, sym, film, l_eibz, l_timeReversalCheck)
      use m_juDFT
      USE m_types_sym
      CLASS(t_kpts), INTENT(inout):: kpts
      TYPE(t_sym), INTENT(IN)     :: sym
      LOGICAL, INTENT(IN)         :: film, l_eibz
      LOGICAL, INTENT(IN)         :: l_timeReversalCheck

      INTEGER :: n,itet,ntet
      call timestart("init_kpts")
      call kpts%find_gamma()
      IF (kpts%nkptf == 0) CALL gen_bz(kpts, sym, l_timeReversalCheck)

      if(l_eibz) then
         kpts%l_set_eibz = .True.
         allocate(kpts%EIBZ(kpts%nkpt))
         !$OMP PARALLEL do default(none) private(n) shared(kpts, sym)
         do n = 1,kpts%nkpt
            call kpts%EIBZ(n)%init(kpts, sym, n)
         enddo
         !$OMP END PARALLEL DO
      end if

      call timestop("init_kpts")
   END SUBROUTINE init_kpts

   subroutine find_gamma_kpts(kpts)
      implicit none 
      class(t_kpts), INTENT(inout):: kpts
      integer :: n 

      kpts%l_gamma = .FALSE.

      DO n = 1, kpts%nkpt
         kpts%l_gamma = kpts%l_gamma .OR. ALL(ABS(kpts%bk(:, n)) < 1E-9)
      ENDDO
   end subroutine find_gamma_kpts

   SUBROUTINE gen_bz(kpts, sym, l_timeReversalCheck)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! gen_bz generates the (whole) Brillouin zone from the          !
      ! (irreducible) k-points given                                  !
      !                                     M.Betzinger (09/07)       !
      !                        Refactored in 2017 by G.M.             !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !     bk     ::    irreducible k-points
      !     nkpt   ::    number of irr. k-points
      !     bkf    ::    all k-points
      !     nkptf  ::    number of all k-points
      !     bkp    ::    k-point parent
      !     bksym  ::    symmetry operation, that connects the parent
      !                  k-point with the current one
      USE m_juDFT
      USE m_types_sym
      TYPE(t_kpts), INTENT(INOUT) :: kpts
      TYPE(t_sym), INTENT(IN)     :: sym
      LOGICAL, INTENT(IN)         :: l_timeReversalCheck
      !  - local scalars -
      INTEGER                  ::  iop, nkptfCheck

      !  - local arrays -
      INTEGER, ALLOCATABLE     ::  iarr(:)
      REAL                     ::  rrot(3, 3, 2*sym%nop)
      REAL, ALLOCATABLE        ::  rarr1(:, :)

      INTEGER:: nsym, ID_mat(3, 3)
      call timestart("gen_bz")

      nsym = sym%nop
      id_mat = 0
      id_mat(1, 1) = 1; id_mat(2, 2) = 1; id_mat(3, 3) = 1
      IF (ANY(sym%mrot(:, :, 1) .NE. id_mat)) CALL judft_error("Identity must be first symmetry operation", calledby="gen_bz")

      IF(l_timeReversalCheck) THEN
         ALLOCATE (kpts%bkf(3, nsym*kpts%nkpt))
         ALLOCATE (kpts%bkp(nsym*kpts%nkpt))
         ALLOCATE (kpts%bksym(nsym*kpts%nkpt))

         ! Generate symmetry operations in reciprocal space
         DO iop = 1, nsym
            rrot(:, :, iop) = TRANSPOSE(sym%mrot(:, :, sym%invtab(iop)))
         END DO

         CALL gen_bz_internal(kpts, sym, rrot, nsym)

         nkptfCheck = kpts%nkptf

         DEALLOCATE (kpts%bksym)
         DEALLOCATE (kpts%bkp)
         DEALLOCATE (kpts%bkf)
      END IF

      IF (.NOT. sym%invs) nsym = 2*sym%nop

      ALLOCATE (kpts%bkf(3, nsym*kpts%nkpt))
      ALLOCATE (kpts%bkp(nsym*kpts%nkpt))
      ALLOCATE (kpts%bksym(nsym*kpts%nkpt))

      ! Generate symmetry operations in reciprocal space
      DO iop = 1, nsym
         IF (iop .LE. sym%nop) THEN
            rrot(:, :, iop) = TRANSPOSE(sym%mrot(:, :, sym%invtab(iop)))
         ELSE
            rrot(:, :, iop) = -rrot(:, :, iop - sym%nop)
         END IF
      END DO

      CALL gen_bz_internal(kpts, sym, rrot, nsym)
      
      IF(l_timeReversalCheck) THEN
         IF(nkptfCheck.NE.kpts%nkptf) THEN
            CALL juDFT_warn("k-point set is not compatible to missing time-reversal symmetry in calculation.",calledby="gen_bz")
         END IF
      END IF
      ! Reallocate bkf, bkp, bksym
      ALLOCATE (iarr(kpts%nkptf))
      iarr = kpts%bkp(:kpts%nkptf)
      DEALLOCATE (kpts%bkp)
      ALLOCATE (kpts%bkp(kpts%nkptf))
      kpts%bkp = iarr
      iarr = kpts%bksym(:kpts%nkptf)
      DEALLOCATE (kpts%bksym)
      ALLOCATE (kpts%bksym(kpts%nkptf))
      kpts%bksym = iarr
      DEALLOCATE (iarr)
      ALLOCATE (rarr1(3, kpts%nkptf))
      rarr1 = kpts%bkf(:, :kpts%nkptf)
      DEALLOCATE (kpts%bkf)
      ALLOCATE (kpts%bkf(3, kpts%nkptf))
      kpts%bkf = rarr1
      DEALLOCATE (rarr1)
      call timestop("gen_bz")
   END SUBROUTINE gen_bz
   
   SUBROUTINE gen_bz_internal(kpts, sym, rrot, nsym)
      USE m_juDFT
      USE m_types_sym
      TYPE(t_kpts), INTENT(INOUT) :: kpts
      TYPE(t_sym), INTENT(IN)     :: sym
      REAL, INTENT(IN)            :: rrot(3, 3, 2*sym%nop)
      INTEGER, INTENT(IN)         :: nsym
      !  - local scalars -
      INTEGER                 ::  ic, iop, ikpt, ikpt1
      INTEGER                 ::  binX, binY, binZ, maxBinSize, iParent

      !  - local arrays -
      INTEGER, ALLOCATABLE     ::  nkptInBin(:,:,:)
      INTEGER, ALLOCATABLE     ::  kptParentBins(:,:,:,:)
      REAL                     ::  rotkpt(3)
      REAL, PARAMETER          ::  eps = 1e-5
      REAL, PARAMETER          ::  sameKPTEps = 1.0e-6

      !Add existing vectors to list of full vectors
      !For a DFPT test calculation, this broke --> set additional .FALSE.
      IF (.FALSE..AND.((kpts%kptsKind.EQ.KPTS_KIND_MESH).OR.(kpts%kptsKind.EQ.KPTS_KIND_SPEX_MESH)).AND.(.NOT.ANY(kpts%nkpt3(:).EQ.0))) THEN
         ALLOCATE (nkptInBin(-(kpts%nkpt3(1)+1):(kpts%nkpt3(1)+1),-(kpts%nkpt3(2)+1):(kpts%nkpt3(2)+1),-(kpts%nkpt3(3)+1):(kpts%nkpt3(3)+1)))
         nkptInBin = 0
         DO ikpt = 1, kpts%nkpt
            binX = FLOOR((kpts%bk(1,ikpt)*kpts%nkpt3(1))+eps)
            binY = FLOOR((kpts%bk(2,ikpt)*kpts%nkpt3(2))+eps)
            binZ = FLOOR((kpts%bk(3,ikpt)*kpts%nkpt3(3))+eps)
            nkptInBin(binX,binY,binZ) = nkptInBin(binX,binY,binZ) + 1
         END DO
         maxBinSize = MAXVAL(nkptInBin) + 2
         DEALLOCATE (nkptInBin) 
         ALLOCATE (kptParentBins(maxBinSize,-(kpts%nkpt3(1)+1):(kpts%nkpt3(1)+1),-(kpts%nkpt3(2)+1):(kpts%nkpt3(2)+1),-(kpts%nkpt3(3)+1):(kpts%nkpt3(3)+1)))
         kptParentBins = 0
         ic = 0
         DO iop = 1, nsym
            DO ikpt = 1, kpts%nkpt
               rotkpt = MATMUL(rrot(:, :, iop), kpts%to_first_bz(kpts%bk(:, ikpt)))
               !transform back into 1st-BZ (Do not use nint to deal properly with inaccuracies)
               rotkpt = kpts%to_first_bz(rotkpt)
               binX = FLOOR((rotkpt(1)*kpts%nkpt3(1))+eps)
               binY = FLOOR((rotkpt(2)*kpts%nkpt3(2))+eps)
               binZ = FLOOR((rotkpt(3)*kpts%nkpt3(3))+eps)
               DO iParent = 1, maxBinSize
                  IF (kptParentBins(iParent,binX,binY,binZ).EQ.0) THEN
                     ic = ic + 1
                     kptParentBins(iParent,binX,binY,binZ) = ic
                     kpts%bkf(:, ic) = rotkpt
                     kpts%bkp(ic) = ikpt
                     kpts%bksym(ic) = iop
                     EXIT
                  ELSE
                     IF (all(abs(kpts%bkf(:,kptParentBins(iParent,binX,binY,binZ)) - rotkpt) < sameKPTEps)) EXIT
                  END IF
               END DO
               IF (iParent.GT.maxBinSize) THEN
                  WRITE(*,*) 'bin size: ', maxBinSize
                  CALL juDFT_error("Bin size too small", calledby='types_kpts%gen_bz')
               END IF
            END DO
         END DO
         DEALLOCATE (kptParentBins)
      ELSE
         ic = 0
         DO iop = 1, nsym
            DO ikpt = 1, kpts%nkpt
               rotkpt = MATMUL(rrot(:, :, iop), kpts%to_first_bz(kpts%bk(:, ikpt)))
               !transform back into 1st-BZ (Do not use nint to deal properly with inaccuracies)
               rotkpt = kpts%to_first_bz(rotkpt)
               DO ikpt1 = 1, ic
                  IF (all(abs(kpts%bkf(:, ikpt1) - rotkpt) < sameKPTEps)) EXIT
               END DO

               IF (ikpt1 > ic) THEN !new point
                  ic = ic + 1
                  kpts%bkf(:, ic) = rotkpt
                  kpts%bkp(ic) = ikpt
                  kpts%bksym(ic) = iop
               END IF
            END DO
         END DO
      END IF
      
      kpts%nkptf = ic
   END SUBROUTINE gen_bz_internal

   function nkpt3_kpts(kpts) result(nkpt3)
      implicit none
      class(t_kpts), intent(in) :: kpts
      integer                   :: nkpt3(3)

      integer :: ikpt
      real    :: k(3)

      nkpt3 = 0

      do ikpt = 1, kpts%nkptf
         k = kpts%bkf(:, ikpt)
#ifdef CPP_EXPLICIT_HYB
         write (*, "(I5,A,F8.4,F8.4,F8.4)") ikpt, ": ", kpts%bkf(:, ikpt)
#endif
         if (abs(k(2)) < 1e-10 .and. abs(k(3)) < 1e-10) nkpt3(1) = nkpt3(1) + 1
         if (abs(k(1)) < 1e-10 .and. abs(k(3)) < 1e-10) nkpt3(2) = nkpt3(2) + 1
         if (abs(k(1)) < 1e-10 .and. abs(k(2)) < 1e-10) nkpt3(3) = nkpt3(3) + 1
      enddo
   end function nkpt3_kpts

   subroutine calc_nkpt_EIBZ(eibz, kpts, sym, nk)
      USE m_types_sym
      USE m_juDFT
      IMPLICIT NONE
      class(t_eibz), intent(inout)  :: eibz
      type(t_kpts),  INTENT(IN)     :: kpts
      TYPE(t_sym),   INTENT(IN)     :: sym
      integer, intent(in)           :: nk

!     - local scalars -
      INTEGER               ::  isym, ic, iop, ikpt, ikpt1
      INTEGER               ::  nsymop, nrkpt
!     - local arrays -
      INTEGER               ::  i
      INTEGER               ::  neqvkpt(kpts%nkptf), list(kpts%nkptf), parent(kpts%nkptf), &
                               symop(kpts%nkptf)
      INTEGER, ALLOCATABLE  ::  psym(:)
      REAL                  ::  rotkpt(3), rrot(3, 3, sym%nsym)
      allocate (psym(sym%nsym))

      ! calculate rotations in reciprocal space
      DO isym = 1, sym%nsym
         IF (isym <= sym%nop) THEN
            rrot(:, :, isym) = transpose(sym%mrot(:, :, sym%invtab(isym)))
         ELSE
            rrot(:, :, isym) = -rrot(:, :, isym - sym%nop)
         END IF
      END DO

      ! determine little group of k., i.e. those symmetry operations
      ! which keep bk(:,nk,nw) invariant
      ! nsymop :: number of such symmetry-operations
      ! psym   :: points to the symmetry-operation

      call calc_psym_nsymop(kpts, sym, nk, psym, nsymop)

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k

      neqvkpt = 0

      !       list = [(ikpt-1, ikpt=1,nkpt) ]
      DO ikpt = 1, kpts%nkptf
         list(ikpt) = ikpt - 1
      END DO

      DO ikpt = 2, kpts%nkptf
         DO iop = 1, nsymop
            
            rotkpt = matmul(rrot(:, :, psym(iop)), kpts%bkf(:, ikpt))

            !transfer rotkpt into BZ
            rotkpt = kpts%to_first_bz(rotkpt)

            !determine number of rotkpt
            nrkpt = 0
            DO ikpt1 = 1, kpts%nkptf
               IF (all(abs(rotkpt - kpts%bkf(:, ikpt1)) <= 1E-06)) THEN
                  nrkpt = ikpt1
                  EXIT
               END IF
            END DO
            IF (nrkpt == 0) call judft_error('symm: Difference vector not found !')
            IF (list(nrkpt) /= 0) THEN
               list(nrkpt) = 0
               neqvkpt(ikpt) = neqvkpt(ikpt) + 1
               parent(nrkpt) = ikpt
               symop(nrkpt) = psym(iop)
            END IF

            if(all(list == 0)) exit
         END DO
      END DO


      ! for the Gamma-point holds:
      parent(1) = 1
      neqvkpt(1) = 1

      ! determine number of members in the EIBZ(k)
      ic = 0
      DO ikpt = 1, kpts%nkptf
         IF (parent(ikpt) == ikpt) ic = ic + 1
      END DO
      EIBZ%nkpt = ic
   END subroutine calc_nkpt_EIBZ

   subroutine calc_pointer_EIBZ(eibz, kpts, sym, nk)
      USE m_types_sym
      implicit none
      class(t_eibz), intent(inout)  :: eibz
      type(t_kpts), INTENT(IN)      :: kpts
      type(t_sym), intent(in)       :: sym
      integer, intent(in)           :: nk

      INTEGER               :: list(kpts%nkptf), parent(kpts%nkptf)
      integer               :: isym, i, nsymop, ic, ikpt, iop, nrkpt
      INTEGER               :: rrot(3, 3, sym%nsym)
      INTEGER, ALLOCATABLE  :: psym(:)
      REAL                  :: rotkpt(3)

      allocate (psym(sym%nsym))
      parent = 0
      DO isym = 1, sym%nsym
         IF (isym <= sym%nop) THEN
            rrot(:, :, isym) = transpose(sym%mrot(:, :, sym%invtab(isym)))
         ELSE
            rrot(:, :, isym) = -rrot(:, :, isym - sym%nop)
         END IF
      END DO

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k

      DO i = 1, kpts%nkptf
         list(i) = i - 1
      END DO

      ! calc numsymop
      call calc_psym_nsymop(kpts, sym, nk, psym, nsymop)

      DO ikpt = 2, kpts%nkptf
         DO iop = 1, nsymop
            rotkpt = matmul(rrot(:, :, psym(iop)), kpts%bkf(:, ikpt))

            !determine number of rotkpt
            nrkpt = kpts%get_nk(rotkpt)
            IF(nrkpt == 0) call judft_error('symm: Difference vector not found !')


            IF(list(nrkpt) /= 0) THEN
               list(nrkpt) = 0
               parent(nrkpt) = ikpt
            END IF
            IF(all(list == 0)) EXIT

         END DO
      END DO

      ! for the Gamma-point holds:
      parent(1) = 1
      IF(ALLOCATED(eibz%pointer)) DEALLOCATE(eibz%pointer)
      allocate(eibz%pointer(eibz%nkpt), source=0)
      ic = 0
      DO ikpt = 1, kpts%nkptf
         IF(parent(ikpt) == ikpt) THEN
            ic = ic + 1
            eibz%pointer(ic) = ikpt
         END IF
      END DO
   end subroutine calc_pointer_EIBZ

   subroutine calc_psym_nsymop(kpts, sym, nk, psym, nsymop)
      USE m_types_sym
      implicit none
      class(t_kpts), intent(in) :: kpts
      type(t_sym), intent(in)   :: sym
      integer, intent(in)       :: nk
      INTEGER, intent(inout)    :: psym(:)
      integer, intent(inout)    :: nsymop

      integer :: ic, iop, isym
      REAL    :: rotkpt(3)
      real    :: rrot(3, 3, sym%nsym)

      DO isym = 1, sym%nsym
         IF (isym <= sym%nop) THEN
            rrot(:, :, isym) = transpose(sym%mrot(:, :, sym%invtab(isym)))
         ELSE
            rrot(:, :, isym) = -rrot(:, :, isym - sym%nop)
         END IF
      END DO

      ic = 0
      DO iop = 1, sym%nsym
         rotkpt = matmul(rrot(:, :, iop), kpts%bkf(:, nk))

         !transfer rotkpt into BZ
         rotkpt = kpts%to_first_bz(rotkpt)

         !check if rotkpt is identical to bk(:,nk)
         IF (maxval(abs(rotkpt - kpts%bkf(:, nk))) <= 1E-07) THEN
            ic = ic + 1
            psym(ic) = iop
         END IF
      END DO
      nsymop = ic
   end subroutine calc_psym_nsymop

   SUBROUTINE calcCommonFractions(kpts,commonFractions)

      USE m_constants

      IMPLICIT NONE

      CLASS(t_kpts), INTENT(IN) :: kpts
      REAL, INTENT(INOUT)    :: commonFractions(3)

      INTEGER, PARAMETER :: upperBound = 50
      INTEGER            :: i, j, ikpt
      REAL               :: temp
      LOGICAL            :: l_CommonFraction

      commonFractions(:) = -1

!      IF(kpts%kptsKind.NE.KPTS_KIND_MESH) THEN
!         commonFractions(:) = -1
!         RETURN
!      END IF

      DO j = 1, 3
         DO i = 2, upperBound
            l_CommonFraction = .TRUE.
            DO ikpt = 1, kpts%nkpt
               temp = i * kpts%bk(j,ikpt)
               IF (ABS(temp - NINT(temp)).GT.1.0e-8) THEN
                  l_CommonFraction = .FALSE.
                  EXIT
               END IF
            END DO
            IF(l_CommonFraction) THEN
               commonFractions(j) = i
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE

END MODULE m_types_kpts
