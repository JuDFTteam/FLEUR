!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_sliceplot
    TYPE t_sliceplot
     LOGICAL :: iplot=.false.
     LOGICAL :: slice=.false.
     LOGICAL :: plpot=.false.
     INTEGER :: kk=0
     INTEGER :: nnne=0
     REAL    :: e1s=0.
     REAL    :: e2s=0.
   contains
     procedure :: read_xml
  END TYPE t_sliceplot

contains
   subroutine read_xml(slicelplot,xml)
     use m_types_xml
     class(t_sliceplot),INTENT(OUT)::sliceplot
     type(t_xml),INTENT(IN)::xml

     character(len=200)::xpatha
     integer::numberNodes

     xPathA = '/fleurInput/output'
     numberNodes = xml%GetNumberOfNodes(xPathA)

     IF (numberNodes.EQ.1) THEN
        siceplot%slice = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@slice'))
     endif
     xPathA = '/fleurInput/output/plotting'
     numberNodes = xml%GetNumberOfNodes(xPathA)

     IF (numberNodes.EQ.1) THEN
        sliceplot%iplot = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@iplot'))
        sliceplot%plpot = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plplot'))
     END IF

     xPathA = '/fleurInput/output/chargeDensitySlicing'
     numberNodes = xmlGetNumberOfNodes(xPathA)

     IF ((sliceplot%slice).AND.(numberNodes.EQ.0)) THEN
        CALL juDFT_error("slice is true but chargeDensitySlicing parameters are not set!", calledby = "r_inpXML")
     END IF

     IF (numberNodes.EQ.1) THEN
        sliceplot%kk = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numkpt'))
        sliceplot%e1s = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minEigenval'))
        sliceplot%e2s = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxEigenval'))
        sliceplot%nnne = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nnne'))
     END IF
   end subroutine read_xml
     
END MODULE m_types_sliceplot
