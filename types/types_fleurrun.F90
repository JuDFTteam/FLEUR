!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_fleurrun
  IMPLICIT NONE
  TYPE t_run
     INTEGER:: itmax
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
   CONTAINS
     PROCEDURE,PASS :: read_xml=>read_xml_run
  END TYPE t_run

CONTAINS
  SUBROUTINE read_xml_run(run)
    USE m_xmlIntWrapFort
    USE m_calculator
    CLASS(t_run),INTENT(INOUT)::run
    
    run%itmax = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@itmax'))
    run%minDistance = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@minDistance'))
    run%swsp = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@swsp'))
    run%lflip = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@lflip'))
    run%skip_pot=evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@skip_pot'))
    run%skip_eigen=evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@skip_eigen'))

    run%score = .FALSE.
    run%secvar = .FALSE.
    run%l_bmt=.FALSE.
    IF (xmlGetNumberOfNodes('/fleurInput/output').EQ.1) THEN
       IF (xmlGetNumberOfNodes('/fleurInput/output/plotting')==1) &
            run%score = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/plotting/@score'))
       IF (xmlGetNumberOfNodes('/fleurInput/output/specialOutput').EQ.1) THEN
          run%l_bmt = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/specialOutput/@bmt'))
          run%eonly = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/specialOutput//@eonly'))
       ENDIF
    END IF
    IF (xmlGetNumberOfNodes('/fleurInput/calculationSetup/expertModes')==1)&
         run%secvar = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/expertModes/@secvar'))

  END SUBROUTINE read_xml_run
END MODULE m_types_fleurrun

  

  
