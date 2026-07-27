!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_libxc_xctyp
   !Parses the inpgen "LibXC: ..." xctyp specification (raw xcpot%inbuild_name)
   !into either name-based (Exch/Cor) or id-based (ExchID/CorID) functional
   !selection, plus an optional auxiliary GGA (AuxExchID/AuxCorID) for the
   !MetaGGA radial basis. Shared between src/fleur/io/w_inpXML.f90 (fleur and
   !inpgen2) and src/tools/inpgen3/w_inpXML.f90 (inpgen3).
   IMPLICIT NONE
CONTAINS

   SUBROUTINE parse_libxc_xctyp(inbuild_name, l_useID, xName, cName, xID, cID, l_aux, auxXID, auxCID)
      USE m_judft
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN)  :: inbuild_name
      LOGICAL,          INTENT(OUT) :: l_useID
      CHARACTER(len=*), INTENT(OUT) :: xName, cName
      INTEGER,          INTENT(OUT) :: xID, cID
      LOGICAL,          INTENT(OUT) :: l_aux
      INTEGER,          INTENT(OUT) :: auxXID, auxCID

      CHARACTER(len=LEN(inbuild_name)) :: body, token, key, value
      INTEGER :: colonIndex, commaIndex, prefixIndex, ios
      LOGICAL :: haveExch, haveCor, haveExchID, haveCorID, haveAuxExchID, haveAuxCorID

      xName = ''; cName = ''
      xID = -1; cID = -1; auxXID = 0; auxCID = 0
      haveExch = .FALSE.; haveCor = .FALSE.
      haveExchID = .FALSE.; haveCorID = .FALSE.
      haveAuxExchID = .FALSE.; haveAuxCorID = .FALSE.

      prefixIndex = INDEX(inbuild_name, ':')
      IF (prefixIndex == 0) CALL judft_error("Malformed LibXC specification: "//TRIM(inbuild_name), &
         calledby="parse_libxc_xctyp", &
         hint='Use e.g. xctyp="LibXC: Exch: LDA_X, Cor: LDA_C_PW" or xctyp="LibXC: ExchID: 645, CorID: 642"')
      body = ADJUSTL(inbuild_name(prefixIndex+1:))

      DO WHILE (LEN_TRIM(body) > 0)
         commaIndex = INDEX(body, ',')
         IF (commaIndex == 0) THEN
            token = body
            body = ''
         ELSE
            token = body(1:commaIndex-1)
            body = ADJUSTL(body(commaIndex+1:))
         ENDIF

         colonIndex = INDEX(token, ':')
         IF (colonIndex == 0) CALL judft_error("Malformed LibXC token: "//TRIM(token), calledby="parse_libxc_xctyp", &
            hint="Every entry must have the form Key: Value, e.g. ExchID: 645")
         key = ADJUSTL(token(1:colonIndex-1))
         value = ADJUSTL(token(colonIndex+1:))
         IF (LEN_TRIM(value) == 0) CALL judft_error("Missing value for key "//TRIM(key)//" in LibXC specification", &
            calledby="parse_libxc_xctyp")

         SELECT CASE (TRIM(key))
         CASE ('Exch')
            xName = TRIM(value); haveExch = .TRUE.
         CASE ('Cor')
            cName = TRIM(value); haveCor = .TRUE.
         CASE ('ExchID')
            READ(value, *, iostat=ios) xID
            IF (ios /= 0) CALL judft_error("Could not read an integer from ExchID: "//TRIM(value), calledby="parse_libxc_xctyp")
            haveExchID = .TRUE.
         CASE ('CorID')
            READ(value, *, iostat=ios) cID
            IF (ios /= 0) CALL judft_error("Could not read an integer from CorID: "//TRIM(value), calledby="parse_libxc_xctyp")
            haveCorID = .TRUE.
         CASE ('AuxExchID')
            READ(value, *, iostat=ios) auxXID
            IF (ios /= 0) CALL judft_error("Could not read an integer from AuxExchID: "//TRIM(value), calledby="parse_libxc_xctyp")
            haveAuxExchID = .TRUE.
         CASE ('AuxCorID')
            READ(value, *, iostat=ios) auxCID
            IF (ios /= 0) CALL judft_error("Could not read an integer from AuxCorID: "//TRIM(value), calledby="parse_libxc_xctyp")
            haveAuxCorID = .TRUE.
         CASE DEFAULT
            CALL judft_error("Unknown key '"//TRIM(key)//"' in LibXC specification", calledby="parse_libxc_xctyp", &
               hint="Allowed keys are Exch, Cor, ExchID, CorID, AuxExchID, AuxCorID")
         END SELECT
      ENDDO

      IF ((haveExch .OR. haveCor) .AND. (haveExchID .OR. haveCorID)) THEN
         CALL judft_error("Specify the LibXC functional either by name (Exch/Cor) or by numeric id (ExchID/CorID), not both", &
            calledby="parse_libxc_xctyp")
      ENDIF
      IF (haveExchID .NEQV. haveCorID) THEN
         CALL judft_error("Both ExchID and CorID must be given together", calledby="parse_libxc_xctyp")
      ENDIF
      IF (haveExch .NEQV. haveCor) THEN
         CALL judft_error("Both Exch and Cor must be given together", calledby="parse_libxc_xctyp")
      ENDIF
      IF (haveAuxExchID .NEQV. haveAuxCorID) THEN
         CALL judft_error("Both AuxExchID and AuxCorID must be given together", calledby="parse_libxc_xctyp")
      ENDIF
      IF (.NOT. (haveExchID .OR. haveExch)) THEN
         CALL judft_error("No exchange/correlation functional given for LibXC, use Exch/Cor or ExchID/CorID", &
            calledby="parse_libxc_xctyp")
      ENDIF

      l_useID = haveExchID
      l_aux = haveAuxExchID

   END SUBROUTINE parse_libxc_xctyp

END MODULE m_libxc_xctyp
