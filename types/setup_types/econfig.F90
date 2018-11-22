!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_econfig
CONTAINS
  SUBROUTINE parse_econfig(inputstring,coreconfig,Occs,Nprnc,kappa,states_provided)
    USE m_judft
    USE m_constants
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)   :: inputstring
    LOGICAL,INTENT(IN)            :: coreconfig
    REAL,INTENT(OUT):: occs(:,:)
    INTEGER,INTENT(OUT):: Nprnc(:),kappa(:)
    INTEGER,INTENT(INOUT):: states_provided

    CHARACTER(:),ALLOCATABLE :: string,token
    INTEGER:: i,j

    string=inputstring
    token = popFirstStringToken(string)
    DO WHILE (token.NE.' ')
       IF (token(1:1).EQ.'[') THEN
          IF (.NOT.coreconfig) CALL judft_error("[] notation only allowed in core configuration",calledby="econfig")
          DO i = 1, 6
             IF (TRIM(ADJUSTL(token)).EQ.nobleGasConfigList_const(i)) THEN
                IF (states_provided+nobleGasNumStatesList_const(i).GT.SIZE(occs,1)) &
                     CALL judft_error('Error: Too many core states provided in xml input file!')
                DO j = states_provided+1, states_provided+nobleGasNumStatesList_const(i)
                   occs(j,:) = coreStateNumElecsList_const(i)
                   Nprnc(j) = coreStateNprncList_const(i)
                   Kappa(j) = coreStateKappaList_const(i)
                END DO
                states_provided = states_provided + nobleGasNumStatesList_const(i)
             END IF
          END DO
       ELSE
          DO i = 1, 29
             IF (TRIM(ADJUSTL(token)).EQ.coreStateList_const(i)) THEN
                states_provided = states_provided + 1
                IF (states_provided.GT.SIZE(occs,1)) &
                     CALL judft_error('Error: Too many core states provided in xml input file!')
                Occs(states_provided,:) = coreStateNumElecsList_const(i)
                Nprnc(states_provided) = coreStateNprncList_const(i)
                Kappa(states_provided) = coreStateKappaList_const(i)
             END IF
          END DO
       END IF
       token = popFirstStringToken(String)
    END DO

  END SUBROUTINE parse_econfig

  SUBROUTINE parse_occupation(xpath,Occs,Nprnc,Kappa,states_provided)
    USE m_judft
    USE m_xmlIntWrapFort
    USE m_calculator
    USE m_constants
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)::xpath
    INTEGER,INTENT(in)::nprnc(:),kappa(:),states_provided
    REAL,INTENT(INOUT)::occs(:,:)

    CHARACTER(len=50):: valueString
    INTEGER          :: nprncTemp,kappaTemp,i
    valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@state')
    nprncTemp = 0
    kappaTemp = 0
    DO i = 1, 29
       IF (TRIM(ADJUSTL(valueString)).EQ.coreStateList_const(i)) THEN
          nprncTemp = coreStateNprncList_const(i)
          kappaTemp = coreStateKappaList_const(i)
       END IF
    END DO
    IF (nprncTemp==0) CALL judft_error("Illegal state definition in "//xpath)
    DO i = 1, states_provided
       IF ((nprncTemp.EQ.Nprnc(i)).AND.(kappaTemp.EQ.Kappa(i))) THEN
          Occs(i,1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@spinUp'))
          Occs(i,2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@spinDown'))
          RETURN
       END IF
    END DO
    CALL judft_error("You can not specify an occupation for a state not explicitly defined")
  END SUBROUTINE parse_occupation

  FUNCTION popFirstStringToken(line) RESULT(firstToken)

    IMPLICIT NONE

    CHARACTER(*), INTENT(INOUT) :: line
    CHARACTER(LEN = LEN(line)) :: firstToken

    INTEGER separatorIndex

    separatorIndex = 0
    line = TRIM(ADJUSTL(line))

    separatorIndex = INDEX(line,' ')
    IF (separatorIndex.LE.1) THEN
       firstToken = ' '
    ELSE
       firstToken = line(1:separatorIndex-1)
       line = line(separatorIndex+1:)
    END IF

  END FUNCTION popFirstStringToken
END MODULE m_econfig
