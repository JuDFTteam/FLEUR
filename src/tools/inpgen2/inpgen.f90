!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

PROGRAM inpgen
!----------------------------------------------------------------------------+
!   Set up a FLEUR inp-file from basic input data; for use and docu please   !
!   refer to inpgen.html (or see http://www.flapw.de/docs/inpgen.html)       !
!                                                                            |
!   The program is based on the input conventions of the FLAIR-code, so that !
!   some compatibility is ensured. The symmetry generator was written by     !
!   M.Weinert and implemented in the FLAIR-code by G.Schneider.              !
!                                                                    gb`02   |
!----------------------------------------------------------------------------+
  USE m_juDFT
  USE m_inpgen_help
  use m_inpgen_version
  use m_fleur_dropxmlschema
  USE m_read_inpgen_input
  USE m_make_crystal
  USE m_make_atomic_defaults
  USE m_make_defaults
  USE m_make_kpoints
  USE m_make_magnetism
  USE m_winpxml
  USE m_xsf_io
  USE m_types_input
  USE m_types_atoms
  USE m_types_cell
  USE m_types_sym
  USE m_types_noco
  USE m_types_vacuum
  USE m_types_banddos
  USE m_types_hybinp
  USE m_types_xcpot_inbuild_nofunction
  USE m_types_forcetheo
  USE m_types_kpts
  USE m_types_gfinp
  USE m_types_hub1inp
  USE m_types_enpara

  USE m_types_sliceplot
  USE m_types_stars
  use m_read_old_inp
  use m_fleurinput_read_xml
  USE m_types_mpinp
  USE m_constants
  USE m_types_xml
  USE m_types_juPhon
  use m_make_sym
  USE m_types_profile

      IMPLICIT NONE

      REAL,    ALLOCATABLE :: atompos(:, :),atomid(:),mag_mom(:,:)
      CHARACTER(len=20), ALLOCATABLE :: atomLabel(:)
      LOGICAL               :: l_fullinput,l_explicit,l_inpxml,l_include(4)

      TYPE(t_input)    :: input
      TYPE(t_atoms)    :: atoms
      TYPE(t_cell)     :: cell
      TYPE(t_sym)      :: sym
      TYPE(t_noco)     :: noco
      TYPE(t_vacuum)   :: vacuum
      TYPE(t_banddos)  :: banddos
      TYPE(t_mpinp)    :: mpinp
      TYPE(t_hybinp)   :: hybinp
      TYPE(t_xcpot_inbuild_nf)::xcpot
      TYPE(t_enpara)   :: enpara
      TYPE(t_forcetheo):: forcetheo
      TYPE(t_kpts), ALLOCATABLE :: kpts(:)

      TYPE(t_sliceplot):: sliceplot
      TYPE(t_stars)    :: stars
      TYPE(t_gfinp)    :: gfinp
      TYPE(t_hub1inp)  :: hub1inp
      TYPE(t_enparaXML):: enparaxml
      TYPE(t_juPhon)   :: juPhon
      TYPE(t_profile)  :: profile

      INTEGER            :: idum, kptsUnit, inpOldUnit, ios, inpgenIUnit
      INTEGER            :: iKpts, numKpts, numKptsPath, numNodes, numAddKptsSets, iPoint
      CHARACTER(len=100) :: filename, filename_add
      CHARACTER(len=200) :: xPath
      CHARACTER(len=800) :: line
      CHARACTER(LEN=40)  :: kptsSelection(3)
      CHARACTER(LEN=200) :: tempString, kptsComment
      CHARACTER(len=40), ALLOCATABLE  :: kpts_str(:)
      CHARACTER(len=40), ALLOCATABLE  :: kptsName(:)
      CHARACTER(len=500), ALLOCATABLE :: kptsPath(:)
      INTEGER, ALLOCATABLE :: kptsBZintegration(:)
      LOGICAL, ALLOCATABLE :: kptsGamma(:)
      LOGICAL, ALLOCATABLE :: l_kptsInitialized(:)
      LOGICAL            :: l_exist, l_addPath, l_check, l_oldinpXML

      TYPE(t_xml)::xml

      INTERFACE
       FUNCTION dropDefaultEConfig() BIND(C, name="dropDefaultEconfig")
         USE iso_c_binding
         INTEGER(c_int) dropDefaultEConfig
       END FUNCTION dropDefaultEConfig

       FUNCTION dropDefault2EConfig() BIND(C, name="dropDefault2EConfig")
         USE iso_c_binding
         INTEGER(c_int) dropDefault2EConfig
       END FUNCTION dropDefault2EConfig

       FUNCTION dropOxidesValEConfig() BIND(C, name="dropOxidesValidationEConfig")
         USE iso_c_binding
         INTEGER(c_int) dropOxidesValEConfig
       END FUNCTION dropOxidesValEConfig

       FUNCTION dropProfiles() BIND(C, name="dropProfiles")
         USE iso_c_binding
         INTEGER(c_int) dropProfiles
       END FUNCTION dropProfiles
      END INTERFACE

      CALL judft_init(oUnit,.FALSE.)

      !Start program and greet user
      CALL inpgen_help()
      call inpgen_version()
      call fleur_dropxmlschema()
      l_explicit=judft_was_argument("-explicit")

      filename_add = ""
      IF (judft_was_argument("-add_name")) filename_add = TRIM(judft_string_for_argument("-add_name"))//"_"

      INQUIRE(file='profile.config',exist=l_exist)
      IF (.NOT.l_exist) idum=dropProfiles()

      OPEN(oUnit,file='out')

      INQUIRE(file=TRIM(filename_add)//'inp.xml',exist=l_inpxml)
      IF (l_inpxml.AND..NOT.(judft_was_argument("-inp.xml").or.judft_was_argument("-overwrite")))&
           CALL judft_error(TRIM(filename_add)//"inp.xml exists and can not be overwritten")

      numKpts = 1
      numKptsPath = 0
      l_addPath = .FALSE.
      IF (judft_was_argument("-inp")) THEN
      ELSEIF (judft_was_argument("-inp.xml")) THEN
         !not yet
         l_fullinput = .TRUE.
         CALL xml%init(filename_add,l_fullinput)
         numKpts = xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/kPointLists/kPointList')
         DO iKpts = 1, numKpts
            xPath = ''
            WRITE (xPath, "(a,i0,a)") '/fleurInput/cell/bzIntegration/kPointLists/kPointList[', iKpts, ']/@type'
            numNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPath)))
            IF(numNodes.EQ.1) THEN
               IF(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))).EQ.'path') numKptsPath = numKptsPath + 1
            END IF
         END DO
         CALL xml%FreeResources()
      ELSEIF(judft_was_argument("-f")) THEN
         !read the input
         CALL peekInpgenInput(numKpts, numKptsPath)
         IF(numKpts.EQ.0) numKpts = 1
         IF(numKptsPath.EQ.0) THEN
            l_addPath = .TRUE.
            numKpts = numKpts + 1
         END IF
      ELSE
         CALL judft_error("You should either specify -inp,-inp.xml or -f command line options. Check -h if unsure")
      ENDIF

      numAddKptsSets = 0
      IF(juDFT_was_argument("-kpt")) THEN
         numAddKptsSets = 1
         numKpts = numKpts + numAddKptsSets
      END IF

      ALLOCATE(kpts(numKpts))
      ALLOCATE(kpts_str(numKpts))
      ALLOCATE(kptsName(numKpts))
      ALLOCATE(kptsPath(numKpts))
      ALLOCATE(kptsGamma(numKpts))
      ALLOCATE(l_kptsInitialized(numKpts))
      ALLOCATE(kptsBZintegration(numKpts))
      kpts_str(:)=""
      kptsPath(:)=""
      kptsName(:)=""
      kptsSelection(:) = ''
      l_kptsInitialized(:) = .TRUE.
      kptsBZintegration = BZINT_METHOD_HIST
      kptsGamma = .FALSE.

      CALL profile%init()
      IF (judft_was_argument("-precise")) THEN
         WRITE(*,*) 'NOTE: right now the "-precise" option is under development and experimental.'
         CALL profile%load("precise")
      ELSE IF (judft_was_argument("-profile")) THEN
         WRITE(*,*) 'NOTE: right now the "-profile" option is under development and experimental.'
         CALL profile%load(TRIM(ADJUSTL(judft_string_for_argument("-profile"))))
      END IF

      IF(profile%atomSetup.EQ."oxides_validation") THEN
         INQUIRE(file='oxides_validation.econfig',exist=l_exist)
         IF (.NOT.l_exist) idum=dropOxidesValEconfig()
      ELSE IF (profile%atomSetup.EQ."default2") THEN
         INQUIRE(file='default2.econfig',exist=l_exist)
         IF (.NOT.l_exist) idum=dropDefault2Econfig()
      ELSE
         INQUIRE(file='default.econfig',exist=l_exist)
         IF (.NOT.l_exist) idum=dropDefaultEconfig()
      END IF

      IF (judft_was_argument("-inp")) THEN
         l_kptsInitialized(:) = .FALSE.
         call read_old_inp(input,atoms,cell,stars,sym,noco,vacuum,forcetheo,&
              sliceplot,banddos,enpara,xcpot,kpts(1),hybinp)
         l_fullinput=.TRUE.
      ELSEIF (judft_was_argument("-inp.xml")) THEN
         !not yet
         l_fullinput=.true. !will be set to false if old inp.xml is read
         l_oldinpXML=.true.
         call Fleurinput_read_xml(0,filename_add,cell,sym,atoms,input,noco,vacuum,sliceplot=Sliceplot,banddos=Banddos,&
                                  hybinp=hybinp, xcpot=Xcpot,kptsSelection=kptsSelection,&
                                  kptsArray=kpts,enparaXML=enparaXML,old_version=l_oldinpXML)
         Call Cell%Init(Dot_product(Atoms%Volmts(:),Atoms%Neq(:)))
         call atoms%init(cell)
         Call Sym%Init(Cell,Input%Film)
         CALL xcpot%init(atoms%ntype)
         CALL enpara%init_enpara(atoms,input%jspins,input%film,enparaXML)
      ELSEIF(judft_was_argument("-f")) THEN
         !read the input
         l_kptsInitialized(:) = .FALSE.
         ALLOCATE (sliceplot%plot(1))
         CALL read_inpgen_input(profile,atompos,atomid,mag_mom,atomlabel,kpts_str,kptsName,kptsPath,kptsBZintegration,&
                                kptsGamma,input,sym,noco,vacuum,stars,xcpot,cell,hybinp)


         IF(input%film) sliceplot%plot(1)%zero(3) = -0.5
         IF (l_addPath) THEN
            l_check = .TRUE.
            CALL add_special_points_default(kpts(numKpts),input%film,cell,l_check)
            IF(l_check) THEN
               kpts_str(numKpts-numAddKptsSets) = 'band=240'
               kptsPath(numKpts-numAddKptsSets) = 'default'
               WRITE(kptsName(numKpts-numAddKptsSets),'(a,i0)') "path-", numKpts-numAddKptsSets
            ELSE
               numKpts = numKpts - 1
               WRITE(*,*) 'No default k-point path for band structures for this unit cell type available:'
               WRITE(*,*) 'None is generated by default. If needed a user specific path has to be specified.'
            END IF
         END IF
         l_fullinput=.FALSE.
      ELSE
         CALL judft_error("You should either specify -inp,-inp.xml or -f command line options. Check -h if unsure")
      ENDIF
      IF (.NOT.l_fullinput) THEN
         !First we determine the spacegoup and map the atoms to groups
         CALL make_crystal(input%film,atomid,atompos,mag_mom,atomlabel,vacuum%dvac,noco,cell,sym,atoms)

         !Generate magnetic settings
         CALL make_magnetism(input,noco,atoms,mag_mom)

         !All atom related parameters are set here. Note that some parameters might
         !have been set in the read_input call before by adding defaults to the atompar module
         CALL make_atomic_defaults(input,vacuum,profile,cell ,atoms,enpara)

         !Set all defaults that have not been specified before or can not be specified in inpgen
         CALL make_defaults(atoms,sym,cell,vacuum,input,stars,xcpot,profile,noco,banddos,mpinp,hybinp)
      ENDIF

      IF (numAddKptsSets.EQ.1) THEN
         kptsName(numKpts) = ''
         kpts_str(numKpts) = ''
         kptsPath(numKpts) = ''
         tempString = ''
         tempString = judft_string_for_argument("-kpt")
         IF (INDEX(tempString,"#")>0) THEN
            kptsName(numKpts)=tempString(:INDEX(tempString,"#")-1)
            tempString = TRIM(ADJUSTL(tempString(INDEX(tempString,"#")+1:)))
         END IF
         kpts_str(numKpts) = TRIM(ADJUSTL(tempString))
         IF(INDEX(kpts_str(numKpts),'band=')==1) THEN
            IF (judft_was_argument("-kptsPath")) THEN
               kptsPath(numKpts) = judft_string_for_argument("-kptsPath")
            ELSE
               kptsPath(numKpts) = 'default'
            END IF
         END IF
         IF(kptsName(numKpts).EQ.'') THEN
            IF(kptsPath(numKpts).EQ.'') THEN
               WRITE(kptsName(numKpts),'(a,i0)') "default-", numKpts
            ELSE
               WRITE(kptsName(numKpts),'(a,i0)') "path-", numKpts
            END IF
         END IF
         l_kptsInitialized(numKpts) = .FALSE.
      END IF

      !
      ! k-points can also be modified here
      !
      DO iKpts = 1, numKpts
         IF (l_kptsInitialized(iKpts)) CYCLE
         CALL make_kpoints(kpts(iKpts),cell,sym,hybinp,input%film,noco%l_ss.or.noco%l_soc,&
                           kptsBZintegration(iKpts),kptsGamma(ikpts),kpts_str(iKpts),kptsName(iKpts),kptsPath(iKpts))
         if(hybinp%l_hybrid .and. kpts(iKpts)%kptsKind == KPTS_KIND_MESH) then
            call timestart("Hybrid setup BZ")
            CALL make_sym(sym,cell,atoms,noco ,input,gfinp)
            call kpts(ikpts)%init(sym, input%film,.true.,.FALSE.)
            call timestop("Hybrid setup BZ")
         endif
      END DO

      IF(ALL(kptsSelection(:).EQ.'')) THEN
         kptsSelection(1) = kpts(1)%kptsName ! This may actually be wrong, but it is a backup solution.
         input%bz_integration = kptsBZintegration(1)
         DO iKpts = numKpts, 1, -1
            IF((kpts(iKpts)%kptsKind.EQ.KPTS_KIND_UNSPECIFIED).OR.(kpts(iKpts)%kptsKind.EQ.KPTS_KIND_MESH)) THEN
               kptsSelection(1) = kpts(iKpts)%kptsName
               input%bz_integration = kptsBZintegration(iKpts)
            END IF
            IF(kpts(iKpts)%kptsKind.EQ.KPTS_KIND_PATH) kptsSelection(2) = kpts(iKpts)%kptsName
         END DO
      END IF
      !
      !Now the IO-section
      !
      call determine_includes(l_include)
      IF (.NOT.l_inpxml.or.judft_was_argument("-overwrite").or.l_oldinpXML) THEN
         !the inp.xml file
         !CALL dump_FleurInputSchema()
         filename=TRIM(filename_add)//"inp.xml"
         if (judft_was_argument("-o")) filename=juDFT_string_for_argument("-o")
         INQUIRE(file=filename,exist=l_exist)
         IF(l_exist) CALL system('mv '//trim(filename)//' '//trim(filename)//'_old')
         CALL w_inpxml(&
              atoms,vacuum,input,stars,sliceplot,forcetheo,banddos, juPhon,&
              cell,sym,xcpot,noco ,mpinp,hybinp,kpts,kptsSelection,enpara,gfinp,&
              hub1inp,l_explicit,l_include,filename,filename_add)
         if (.not.l_include(2)) CALL sym%print_XML(99,TRIM(filename_add)//"sym.xml")
      ENDIF

      inpOldUnit = 39
      IF (.NOT.l_include(1)) THEN
         kptsUnit = 38
         INQUIRE(file='kpts.xml',exist=l_exist)
         IF((.NOT.l_exist).AND.judft_was_argument("-inp.xml")) THEN
            CALL system('mv '//TRIM(filename_add)//'inp.xml '//TRIM(filename_add)//'inp_old.xml')
            OPEN (inpOldUnit, file=TRIM(filename_add)//"inp_old.xml", action="read")
            OPEN (kptsUnit, file=TRIM(filename_add)//"inp.xml", action="write")
            ios = 0
            DO WHILE(ios==0)
               READ(inpOldUnit,'(a)',iostat=ios) line
               IF (TRIM(ADJUSTL(line)).EQ.'<kPointLists>') EXIT
               WRITE(kptsUnit,'(a)') TRIM(line)
            END DO
         ELSE
            OPEN (kptsUnit, file=TRIM(filename_add)//"kpts.xml", action="write")
         END IF

         WRITE (kptsUnit, '(a)') "         <kPointLists>"
         DO iKpts = 1, numKpts
            call kpts(iKpts)%find_gamma()
            CALL kpts(iKpts)%print_XML(kptsUnit)
         END DO
         WRITE (kptsUnit, '(a)') "         </kPointLists>"

         IF((.NOT.l_exist).AND.judft_was_argument("-inp.xml")) THEN
            DO WHILE(ios==0)
               READ(inpOldUnit,'(a)',iostat=ios) line
               IF (TRIM(ADJUSTL(line)).EQ.'</kPointLists>') EXIT
            END DO
            DO WHILE(ios==0)
               READ(inpOldUnit,'(a)',iostat=ios) line
               IF(ios.NE.0) EXIT
               WRITE(kptsUnit,'(a)') TRIM(line)
            END DO
         END IF

         CLOSE (kptsUnit)
      END IF

      inpgenIUnit = 57
      IF(judft_was_argument("-f").AND..NOT.juDFT_was_argument("-noInpgenComment")) THEN
         filename = juDFT_string_for_argument("-f")
         OPEN (inpgenIUnit,file=TRIM(filename),action="read")
         OPEN (inpOldUnit, file=TRIM(filename_add)//"inp.xml", action="write", status='old', access='append')
         WRITE(inpOldUnit,'(a)') ''
         WRITE(inpOldUnit,'(a)') '<!--'
         WRITE(inpOldUnit,'(a)') 'Command line when calling inpgen (only for documentation purposes):'
         CALL GET_COMMAND(line)
         WRITE(inpOldUnit,'(a)') TRIM(line)
         WRITE(inpOldUnit,'(a)') ''
         WRITE(inpOldUnit,'(a)') 'Initial (original) inpgen input (only for documentation purposes):'
         ios = 0
         DO WHILE(ios==0)
            READ(inpgenIUnit,'(a)',iostat=ios) line
            IF (ios.EQ.0) WRITE(inpOldUnit,'(a)') TRIM(line)
         END DO
         WRITE(inpOldUnit,'(a)') '-->'
         CLOSE (inpOldUnit)
         CLOSE (inpgenIUnit)
         line = ""
      END IF

100   FORMAT (a20,a15,i10,3x,a)
      WRITE(*,*) 'Stored k-point lists:'
      WRITE(*,*) ''
      WRITE(*,'(a20,a15,a10,3x,a)') 'NAME', 'TYPE', 'NKPT', 'COMMENT'
      WRITE(*,*) '================================================================================'
      DO iKpts = 1, numKpts
         SELECT CASE(kpts(iKpts)%kptsKind)
            CASE (KPTS_KIND_UNSPECIFIED)
               WRITE(*,100) TRIM(ADJUSTL(kpts(iKpts)%kptsName)), 'UNSPECIFIED', kpts(ikpts)%nkpt, ''
            CASE (KPTS_KIND_MESH)
               kptsComment = ''
               WRITE(kptsComment,'(i0,a,i0,a,i0)') kpts(iKpts)%nkpt3(1), ' x ', kpts(iKpts)%nkpt3(2), ' x ', kpts(iKpts)%nkpt3(3)
               WRITE(*,100) TRIM(ADJUSTL(kpts(iKpts)%kptsName)), 'MESH', kpts(iKpts)%nkpt, TRIM(ADJUSTL(kptsComment))
            CASE (KPTS_KIND_PATH)
               kptsComment = ''
               kptsComment = TRIM(ADJUSTL(kpts(iKpts)%specialPointNames(1)))
               DO iPoint = 2, kpts(iKpts)%numSpecialPoints
                  kptsComment = TRIM(ADJUSTL(kptsComment))//' - '//TRIM(ADJUSTL(kpts(iKpts)%specialPointNames(iPoint)))
               END DO
               WRITE(*,100) TRIM(ADJUSTL(kpts(iKpts)%kptsName)), 'PATH', kpts(iKpts)%nkpt, TRIM(ADJUSTL(kptsComment))
            CASE (KPTS_KIND_TRIA_BULK)
               WRITE(*,100) TRIM(ADJUSTL(kpts(iKpts)%kptsName)), 'TRIA-BULK', kpts(ikpts)%nkpt, ''
            CASE (KPTS_KIND_TRIA)
               WRITE(*,100) TRIM(ADJUSTL(kpts(iKpts)%kptsName)), 'TRIA', kpts(ikpts)%nkpt, ''
            CASE (KPTS_KIND_SPEX_MESH)
               kptsComment = ''
               WRITE(kptsComment,'(i0,a,i0,a,i0)') kpts(iKpts)%nkpt3(1), ' x ', kpts(iKpts)%nkpt3(2), ' x ', kpts(iKpts)%nkpt3(3)
               WRITE(*,100) TRIM(ADJUSTL(kpts(iKpts)%kptsName)), 'SPEX-MESH', kpts(iKpts)%nkpt, TRIM(ADJUSTL(kptsComment))
         END SELECT
      END DO
      WRITE(*,*) '================================================================================'

      ! Structure in  xsf-format
      OPEN (55,file="struct.xsf")
      CALL xsf_WRITE_atoms(55,atoms,input%film,cell%amat)
      CLOSE (55)
      CLOSE(oUnit)

      CALL juDFT_end("All done")

    CONTAINS
      SUBROUTINE determine_includes(l_include)
        LOGICAL,INTENT(out)::l_include(4)  !kpts,operations,species,position

        CHARACTER(len=100)::str=''
        LOGICAL           ::incl

        l_include=[.FALSE.,.FALSE.,.TRUE.,.TRUE.]

        IF (judft_was_argument("-inc")) THEN
           str=judft_string_for_argument("-inc")

           DO WHILE(LEN_TRIM(str)>0)
              IF (str(1:1)=='-') THEN
                 incl=.FALSE.
                 str=str(2:)
              ELSE
                 incl=.TRUE.
                 IF (str(1:1)=='+') str=str(2:)
              ENDIF
              SELECT CASE(str(1:1))
              CASE ('k','K')
                 l_include(1)=incl
              CASE ('o','O')
                 l_include(2)=incl
              CASE ('s','S')
                 l_include(3)=incl
              CASE ('p','P')
                 l_include(4)=incl
              CASE ('a','A')
                 l_include(:)=incl
              END SELECT
              IF (INDEX(str,"'")>0) THEN
                 str=str(INDEX(str,"'")+1:)
              ELSE
                 str=""
              END IF
           END DO
        ENDIF

      END SUBROUTINE determine_includes

    END PROGRAM inpgen
