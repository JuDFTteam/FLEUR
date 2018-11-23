!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_job
  IMPLICIT NONE
  TYPE t_job
     INTEGER:: itmax
     LOGICAL:: l_opti
     LOGICAL:: strho
     REAL   :: minDistance
     !optional stuff
     LOGICAL :: swsp
     LOGICAL :: lflip
     LOGICAL :: score
     LOGICAL :: l_bmt
     !Switches controlling the flow of the calculation
     LOGICAL :: secvar
     LOGICAL :: eonly
     LOGICAL :: skip_pot
     LOGICAL :: skip_eigen
     LOGICAL :: l_gw
   CONTAINS
     PROCEDURE,PASS :: read_xml=>read_xml_job
  END TYPE t_job

CONTAINS
  SUBROUTINE read_xml_job(job)
    USE m_xmlIntWrapFort
    USE m_calculator
    CLASS(t_job),INTENT(INOUT)::job
    
    job%itmax = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@itmax'))
    job%minDistance = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@minDistance'))
    job%swsp = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@swsp'))
    job%lflip = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@lflip'))
    job%skip_pot=evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@skip_pot'))
    job%skip_eigen=evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@skip_eigen'))
    job%l_gw=evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@gw_io'))

    job%score = .FALSE.
    job%secvar = .FALSE.
    job%l_bmt=.FALSE.
    IF (xmlGetNumberOfNodes('/fleurInput/output').EQ.1) THEN
       IF (xmlGetNumberOfNodes('/fleurInput/output/plotting')==1) &
            job%score = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/plotting/@score'))
       IF (xmlGetNumberOfNodes('/fleurInput/output/specialOutput').EQ.1) THEN
          job%l_bmt = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/specialOutput/@bmt'))
          job%eonly = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/specialOutput//@eonly'))
       ENDIF
    END IF
    IF (xmlGetNumberOfNodes('/fleurInput/calculationSetup/expertModes')==1)&
         job%secvar = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/expertModes/@secvar'))

  END SUBROUTINE read_xml_job
END MODULE m_types_job

  

  
