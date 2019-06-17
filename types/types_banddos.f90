!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_banddos

   TYPE t_banddos
      LOGICAL :: dos =.false. 
      LOGICAL :: band =.false.
      LOGICAL :: l_mcd =.false.
      LOGICAL :: l_orb =.false.
      LOGICAL :: vacdos =.false.
      INTEGER :: ndir =0
      INTEGER :: orbCompAtom=0
      REAL    :: e1_dos=0.5
      REAL    :: e2_dos=-0.5
      REAL    :: sig_dos=0.015
      REAL    :: e_mcd_lo =-10.0
      REAL    :: e_mcd_up= 0.0
      LOGICAL :: unfoldband =.false.
      INTEGER :: s_cell_x
      INTEGER :: s_cell_y
      INTEGER :: s_cell_z
      REAL    :: alpha,beta,gamma !For orbital decomp. (was orbcomprot)
    contains
      procedure :: read_xml
   END TYPE t_banddos

   subroutine read_xml(banddos,xml)
     use m_types_xml
     class(t_banddos),INTENT(OUT)::banddos
     type(t_xml),INTENT(IN)::xml

     integer::numberNodes

     banddos%band = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@band'))
     
     banddos%dos = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dos'))
     banddos%vacdos = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vacdos'))
     banddos%l_mcd = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mcd'))

     numberNodes = xmlGetNumberOfNodes('/fleurInput/output/densityOfStates')

     IF ((banddos%dos).AND.(numberNodes.EQ.0)) THEN
        CALL juDFT_error("dos is true but densityOfStates parameters are not set!", calledby = "r_inpXML")
     END IF
     
     IF (numberNodes.EQ.1) THEN
        banddos%ndir = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@ndir'))
        banddos%e2_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@minEnergy'))
        banddos%e1_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@maxEnergy'))
        banddos%sig_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/densityOfStates/@sigma'))
     END IF
     IF (banddos%band) THEN
        banddos%dos=.TRUE.
        banddos%ndir = -4
        WRITE(*,*) 'band="T" --> Overriding "dos" and "ndir"!'
     ENDIF

     ! Read in optional magnetic circular dichroism parameters
     numberNodes = xml%GetNumberOfNodes('/fleurInput/output/magneticCircularDichroism')
     
     IF ((banddos%l_mcd).AND.(numberNodes.EQ.0)) THEN
        CALL juDFT_error("mcd is true but magneticCircularDichroism parameters are not set!", calledby = "r_inpXML")
     END IF

     IF (numberNodes.EQ.1) THEN
        banddos%e_mcd_lo = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/magneticCircularDichroism/@energyLo'))
        banddos%e_mcd_up = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/magneticCircularDichroism/@energyUp'))
     END IF

     ! Read in optional parameter for unfolding bandstructure of supercell
     numberNodes = xml%GetNumberOfNodes('/fleurInput/output/unfoldingBand')
     
     IF (numberNodes.EQ.1) THEN
        banddos%unfoldband = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@unfoldband'))
        banddos%s_cell_x = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellX'))
        banddos%s_cell_y = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellY'))
        banddos%s_cell_z = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellZ'))
     END IF
   end subroutine read_xml
     
 END MODULE m_types_banddos
  
