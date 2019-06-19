MODULE m_fleurinput_read_xml
  USE m_types_fleurinput
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleurinput_read_xml(cell,sym,atoms,input,noco,vacuum,field,&
       sliceplot,banddos,hybrid,oneD,coreSpecInput,wann,&
       xcpot,forcetheo_data,kpts,enparaXML)
    USE m_types_xml
    
    TYPE(t_cell),INTENT(OUT)::cell
    TYPE(t_sym),INTENT(OUT)::sym
    TYPE(t_atoms),INTENT(OUT)::atoms
    TYPE(t_input),INTENT(OUT)::input
    TYPE(t_noco),INTENT(OUT)::noco
    TYPE(t_vacuum),INTENT(OUT)::vacuum
    TYPE(t_field),INTENT(OUT)::field
    TYPE(t_sliceplot),INTENT(OUT)::sliceplot
    TYPE(t_banddos),INTENT(OUT)::banddos
    TYPE(t_hybrid),INTENT(OUT)::hybrid
    TYPE(t_oneD),INTENT(OUT)::oneD
    TYPE(t_coreSpecInput),INTENT(OUT)::coreSpecInput
    TYPE(t_wann),INTENT(OUT)::wann
    CLASS(t_xcpot),INTENT(OUT)::xcpot
    TYPE(t_forcetheo_data),INTENT(OUT)::forcetheo_data
    TYPE(t_enparaXML),INTENT(OUT)::enparaXML
    TYPE(t_kpts),INTENT(OUT)::kpts
    
    TYPE(t_xml)::xml

    !Call to init of xml type initialized XML reading and connects to inp.xml
    call xml%init()
    
    !Now read from inp.xml for all datatypes
    CALL cell%read_xml(xml)
    CALL sym%read_xml(xml)
    CALL atoms%read_xml(xml)
    CALL input%read_xml(xml)
    CALL noco%read_xml(xml)
    CALL vacuum%read_xml(xml)
    CALL field%read_xml(xml)
    CALL sliceplot%read_xml(xml)
    CALL banddos%read_xml(xml)
    CALL hybrid%read_xml(xml)
    CALL oneD%read_xml(xml)
    CALL coreSpecInput%read_xml(xml)
    CALL wann%read_xml(xml)
    CALL xcpot%read_xml(xml)
    CALL forcetheo_data%read_xml(xml)
    CALL enparaXML%read_xml(xml)
    CALL kpts%read_xml(xml)

  END SUBROUTINE fleurinput_read_xml
END MODULE m_fleurinput_read_xml
