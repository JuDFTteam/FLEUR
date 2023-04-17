MODULE m_fleurinput_read_xml
  USE m_types_fleurinput
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleurinput_read_xml(xmlOUTFileID,filename_add,cell,sym,atoms,input,noco,vacuum,field,&
       sliceplot,banddos,mpinp,hybinp ,coreSpecInput,wann,&
       xcpot,forcetheo_data,kpts,kptsSelection,kptsArray,enparaXML,gfinp,hub1inp,juPhon,old_version)
    USE m_types_xml
    integer,INTENT(IN)             :: xmlOUTFileID
    CHARACTER(len=100), INTENT(IN) :: filename_add
    TYPE(t_cell),INTENT(OUT),OPTIONAL::cell
    TYPE(t_sym),INTENT(OUT),OPTIONAL::sym
    TYPE(t_atoms),INTENT(OUT),OPTIONAL::atoms
    TYPE(t_input),INTENT(OUT),OPTIONAL::input
    TYPE(t_noco),INTENT(OUT),OPTIONAL::noco
    TYPE(t_vacuum),INTENT(OUT),OPTIONAL::vacuum
    TYPE(t_field),INTENT(OUT),OPTIONAL::field
    TYPE(t_sliceplot),INTENT(OUT),OPTIONAL::sliceplot
    TYPE(t_banddos),INTENT(OUT),OPTIONAL::banddos
    TYPE(t_mpinp), INTENT(OUT), OPTIONAL :: mpinp
    TYPE(t_hybinp),INTENT(OUT),OPTIONAL::hybinp

    TYPE(t_coreSpecInput),INTENT(OUT),OPTIONAL::coreSpecInput
    TYPE(t_wann),INTENT(OUT),OPTIONAL::wann
    CLASS(t_xcpot),INTENT(OUT),OPTIONAL::xcpot
    TYPE(t_forcetheo_data),INTENT(OUT),OPTIONAL::forcetheo_data
    TYPE(t_enparaXML),INTENT(OUT),OPTIONAL::enparaXML
    TYPE(t_kpts),INTENT(OUT),OPTIONAL::kpts
    TYPE(t_kpts),ALLOCATABLE,INTENT(INOUT),OPTIONAL::kptsArray(:)
    TYPE(t_gfinp),INTENT(OUT),OPTIONAL::gfinp
    TYPE(t_hub1inp),INTENT(OUT),OPTIONAL::hub1inp
    TYPE(t_juPhon),INTENT(OUT),OPTIONAL::juPhon
    CHARACTER(LEN=40),INTENT(OUT),OPTIONAL::kptsSelection(3)
    LOGICAL,INTENT(INOUT),OPTIONAL :: old_version

    TYPE(t_xml)::xml

    INTEGER :: numNodes, iNode
    CHARACTER(LEN=40) :: listName, altPurpose
    CHARACTER(LEN=200) :: xPath

    !Call to init of xml type initialized XML reading and connects to inp.xml
    call xml%init(filename_add,old_version)

    !Now read from inp.xml for all datatypes
    if (present(cell)) call cell%read_xml(xml)
    if (present(sym)) call sym%read_xml(xml)
    if (present(atoms)) call atoms%read_xml(xml)
    if (present(input)) call input%read_xml(xml)
    if (present(noco)) call noco%read_xml(xml)
    if (present(vacuum)) call vacuum%read_xml(xml)
    if (present(field)) call field%read_xml(xml)
    if (present(sliceplot)) call sliceplot%read_xml(xml)
    if (present(banddos)) call banddos%read_xml(xml)
    if (present(mpinp)) call mpinp%read_xml(xml)
    if (present(hybinp)) call hybinp%read_xml(xml)
    if (present(coreSpecInput)) call coreSpecInput%read_xml(xml)
    if (present(wann)) call wann%read_xml(xml)
    if (present(xcpot)) call xcpot%read_xml(xml)
    if (present(forcetheo_data)) call forcetheo_data%read_xml(xml)
    if (present(enparaXML)) call enparaXML%read_xml(xml)
    if (present(kpts)) CALL kpts%read_xml(xml)
    if (present(gfinp)) CALL gfinp%read_xml(xml)
    if (present(hub1inp)) CALL hub1inp%read_xml(xml)
    if (present(juPhon)) THEN
        CALL juPhon%read_xml(xml)
        IF (juPhon%kgqmax.LE.0.0) THEN
            juPhon%kgqmax = input%rkmax
        END IF

        IF (juPhon%gqmax.LE.0.0) THEN
            juPhon%gqmax = input%gmax
        END IF
    end if
    IF (present(kptsSelection).and.xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/kPointListSelection')>0) THEN
       kptsSelection(:) = ''
       kptsSelection(1) = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/cell/bzIntegration/kPointListSelection/@listName')))
       numNodes = xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/altKPointList')
       DO iNode = 1 , numNodes
          WRITE (xPath, "(a,i0,a)") '/fleurInput/cell/bzIntegration/altKPointList[',iNode,']'
          altPurpose = ''
          altPurpose = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@purpose')))
          IF (TRIM(ADJUSTL(altPurpose)).EQ.'bands') THEN
             kptsSelection(2) = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@listName')))
          END IF
          IF (TRIM(ADJUSTL(altPurpose)).EQ.'GW') THEN
             kptsSelection(3) = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@listName')))
          END IF
       END DO
    END IF
    IF (PRESENT(kptsArray).and.xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/kPointListSelection')>0) THEN
       numNodes = xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/kPointLists/kPointList')
       IF(.NOT.ALLOCATED(kptsArray)) THEN
          ALLOCATE(kptsArray(numNodes))
       END IF
       DO iNode = 1, numNodes
          CALL kptsArray(iNode)%read_xml_kptsByIndex(xml,iNode)
       END DO
    END IF

    call xml%writexml(xmlOUTFileID)
    call xml%FreeResources()
  END SUBROUTINE fleurinput_read_xml
END MODULE m_fleurinput_read_xml
