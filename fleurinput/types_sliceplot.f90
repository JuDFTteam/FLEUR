!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_sliceplot
  USE m_judft 
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_sliceplot
  TYPE,EXTENDS(t_fleurinput_base) ::t_sliceplot
     LOGICAL :: iplot=.FALSE.
     LOGICAL :: slice=.FALSE.
     LOGICAL :: plpot=.FALSE.
     INTEGER :: kk=0
     INTEGER :: nnne=0
     REAL    :: e1s=0.
     REAL    :: e2s=0.
   CONTAINS
     PROCEDURE :: read_xml=>read_xml_sliceplot
  END TYPE t_sliceplot

CONTAINS
  SUBROUTINE read_xml_sliceplot(this,xml)
    USE m_types_xml
    CLASS(t_sliceplot),INTENT(inOUT)::this
    TYPE(t_xml),INTENT(IN)::xml

    CHARACTER(len=200)::xpatha
    INTEGER::numberNodes

    xPathA = '/fleurInput/output'
    numberNodes = xml%GetNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
       this%slice = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@slice'))
    ENDIF
    xPathA = '/fleurInput/output/plotting'
    numberNodes = xml%GetNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
       this%iplot = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@iplot'))
       this%plpot = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plplot'))
    END IF

    xPathA = '/fleurInput/output/chargeDensitySlicing'
    numberNodes = xml%GetNumberOfNodes(xPathA)

    IF ((this%slice).AND.(numberNodes.EQ.0)) THEN
       CALL juDFT_error("slice is true but chargeDensitySlicing parameters are not set!")
    END IF

    IF (numberNodes.EQ.1) THEN
       this%kk = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numkpt'))
       this%e1s = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minEigenval'))
       this%e2s = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxEigenval'))
       this%nnne = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nnne'))
    END IF
  END SUBROUTINE read_xml_sliceplot

END MODULE m_types_sliceplot
