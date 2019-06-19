!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_coreSpecInput
  USE m_judft
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_coreSpecInput
  ! type for the input to the calculation of the core spectrum (EELS)
  TYPE,EXTENDS(t_fleurinput_base):: t_coreSpecInput
     INTEGER :: verb=0  ! output verbosity
     INTEGER :: atomType  ! atomic type used for calculation of core spectra
     CHARACTER(LEN=1) :: edge  ! edge character (K,L,M,N,O,P)
     INTEGER :: edgeidx(11)  ! l-j edges
     INTEGER :: lx  ! maximum lmax considered in spectra calculation
     REAL :: ek0  ! kinetic energy of incoming electrons
     REAL :: emn  ! energy spectrum lower bound
     REAL :: emx  ! energy spectrum upper bound
     REAL :: ein  ! energy spectrum increment
     INTEGER :: nqphi ! no. of angle-sectors for integral over q vectors
     INTEGER :: nqr   ! no. of radial-sectors for integral over q vectors
     REAL :: alpha_ex  ! maximal angle of incoming electrons
     REAL :: beta_ex   ! maximal (measured) angle of outcoming electrons
     REAL :: I0        ! incoming intensity
   CONTAINS
     PROCEDURE read_xml=>read_xml_corespecinput
  END TYPE t_coreSpecInput
CONTAINS
  SUBROUTINE read_xml_corespecinput(This,xml)
    USE m_types_xml
    CLASS(t_coreSpecInput),INTENT(OUT)::this
    TYPE(t_xml),INTENT(IN)::xml
    
    
    INTEGER:: numberNodes,tempInt,numTokens,i
    LOGICAL::tempBool
    CHARACTER(len=100) :: xPathA,xPathB,valueString
    ! Read in optional core spectrum (EELS) input parameters
    
    xPathA = '/fleurInput/output/coreSpectrum'
    numberNodes = xml%GetNumberOfNodes(xPathA)
    
    IF (numberNodes.EQ.1) THEN
       tempBool = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@verbose'))
       IF(tempBool) this%verb = 1
       this%ek0 = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eKin'))
       this%atomType = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@atomType'))
       this%lx = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lmax'))
       this%edge = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@edgeType')))
       this%emn = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eMin'))
       this%emx = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eMax'))
       tempInt = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numPoints'))
       this%ein = (this%emx - this%emn) / (tempInt - 1.0)
       this%nqphi = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nqphi'))
       this%nqr = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nqr'))
       this%alpha_ex = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@alpha_Ex'))
       this%beta_ex = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@beta_Ex'))
       this%I0 = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@I_initial'))
       xPathB = TRIM(ADJUSTL(xPathA))//'/edgeIndices'
       xPathB = TRIM(ADJUSTL(xPathB))//'/text()'
       valueString = xml%GetAttributeValue(TRIM(ADJUSTL(xPathB)))
       numTokens = xml%countStringTokens(valueString)
       this%edgeidx(:) = 0
       IF(numTokens.GT.SIZE(this%edgeidx)) THEN
          CALL juDFT_error('More EELS edge indices provided than allowed.')
       END IF
       DO i = 1, MAX(numTokens,SIZE(this%edgeidx))
          this%edgeidx(i) = evaluateFirstIntOnly(xml%popFirstStringToken(valueString))
       END DO
    END IF
    
  END SUBROUTINE read_xml_corespecinput
END MODULE m_types_coreSpecInput
