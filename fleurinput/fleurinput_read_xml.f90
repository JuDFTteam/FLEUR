MODULE m_fleurinput_read_xml
  USE m_types_fleurinput
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleurinput_read_xml(cell,sym,atoms,input,noco,vacuum,field,&
       sliceplot,banddos,hybinp,oneD,coreSpecInput,wann,&
       xcpot,forcetheo_data,kpts,enparaXML)
    USE m_types_xml

    TYPE(t_cell),INTENT(OUT),OPTIONAL::cell
    TYPE(t_sym),INTENT(OUT),OPTIONAL::sym
    TYPE(t_atoms),INTENT(OUT),OPTIONAL::atoms
    TYPE(t_input),INTENT(OUT),OPTIONAL::input
    TYPE(t_noco),INTENT(OUT),OPTIONAL::noco
    TYPE(t_vacuum),INTENT(OUT),OPTIONAL::vacuum
    TYPE(t_field),INTENT(OUT),OPTIONAL::field
    TYPE(t_sliceplot),INTENT(OUT),OPTIONAL::sliceplot
    TYPE(t_banddos),INTENT(OUT),OPTIONAL::banddos
    TYPE(t_hybinp),INTENT(OUT),OPTIONAL::hybinp
    TYPE(t_oneD),INTENT(OUT),OPTIONAL::oneD
    TYPE(t_coreSpecInput),INTENT(OUT),OPTIONAL::coreSpecInput
    TYPE(t_wann),INTENT(OUT),OPTIONAL::wann
    CLASS(t_xcpot),INTENT(OUT),OPTIONAL::xcpot
    TYPE(t_forcetheo_data),INTENT(OUT),OPTIONAL::forcetheo_data
    TYPE(t_enparaXML),INTENT(OUT),OPTIONAL::enparaXML
    TYPE(t_kpts),INTENT(OUT),OPTIONAL::kpts

    TYPE(t_xml)::xml

    !Call to init of xml type initialized XML reading and connects to inp.xml
    call xml%init()

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
    if (present(hybinp)) call hybinp%read_xml(xml)
    if (present(oneD)) call oneD%read_xml(xml)
    if (present(coreSpecInput)) call coreSpecInput%read_xml(xml)
    if (present(wann)) call wann%read_xml(xml)
    if (present(xcpot)) call xcpot%read_xml(xml)
    if (present(forcetheo_data)) call forcetheo_data%read_xml(xml)
    if (present(enparaXML)) call enparaXML%read_xml(xml)
    if (present(kpts)) CALL kpts%read_xml(xml)

    call xml%FreeResources()
  END SUBROUTINE fleurinput_read_xml
END MODULE m_fleurinput_read_xml
