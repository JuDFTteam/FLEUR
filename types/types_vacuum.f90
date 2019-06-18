!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_vacuum
  use m_juDFT
  IMPLICIT NONE
   TYPE t_vacuum
      !Stuff for the vacuum
      INTEGER ::nmz=250
      INTEGER ::nmzd=250
      INTEGER ::nmzxy=100
      INTEGER ::nmzxyd=100
      INTEGER :: layerd=1
      INTEGER :: layers=0
      INTEGER :: nvac=2
      INTEGER :: nvacd=2
      REAL :: delz=0.1
      REAL :: dvac=0.0
      INTEGER::nstars=0
      INTEGER:: nstm=0
      REAL :: tworkf=0.0
      REAL :: locx(2)=[0.,0.]
      REAL :: locy(2)=[0.,0.]
      LOGICAL ::starcoeff=.false.
      INTEGER, ALLOCATABLE :: izlay(:, :)
    CONTAINS
      PROCEDURE :: read_xml
   END TYPE t_vacuum
 CONTAINS
   SUBROUTINE read_xml(vacuum,xml)
     USE m_types_xml
     CLASS(t_vacuum),INTENT(OUT)::Vacuum
     TYPE(t_xml),INTENT(IN)::xml
     CHARACTER(len=100)::xpatha
     

     IF (xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1) THEN
     vacuum%dvac = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/cell/filmLattice/@dVac'))
     !vacuum%dtild = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/cell/filmLattice/@dTilda'))

     xPathA = '/fleurInput/output/vacuumDOS'
     IF (xml%GetNumberOfNodes(xpathA)==1) THEN
        vacuum%layers = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@layers'))
        vacuum%starcoeff = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@star'))
        vacuum%nstars = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstars'))
        vacuum%locx(1) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx1'))
        vacuum%locx(2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx2'))
        vacuum%locy(1) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy1'))
        vacuum%locy(2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy2'))
        vacuum%nstm = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstm'))
        vacuum%tworkf = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@tworkf'))
        PRINT *,'Reading of layers not implemented '
     END IF
     vacuum%layerd = vacuum%layers
     ALLOCATE(vacuum%izlay(vacuum%layerd,2))
  ENDIF
   END SUBROUTINE read_xml
 END MODULE m_types_vacuum
