MODULE m_broyd_io

!TODO cdn_IO and getIOMODE IS COMMENTED OUT WHICH IS NOT THE ORIGINAL VERSION!!!!!
USE m_types_hybdat
!USE m_cdn_io

IMPLICIT NONE

CONTAINS

SUBROUTINE readLastIterInAndDiffDen(hybdat,vecLen,nextIter,alpha,inDenVec,diffDenVec, idir, iqpt, iDatom)

   TYPE(t_hybdat), INTENT(IN)  :: hybdat
   INTEGER,        INTENT(IN)  :: vecLen
   integer,        intent(in)  :: idir
   integer,        intent(in)  :: iqpt
   integer,        intent(in)  :: iDatom
   INTEGER,        INTENT(OUT) :: nextIter
   REAL,           INTENT(OUT) :: alpha
   REAL,           INTENT(OUT) :: inDenVec(vecLen), diffDenVec(vecLen)
   character(len=21)           :: filename

   INTEGER             :: mode

   !CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent
   !TODO broyd
   write(filename, '(a7,i1,a4,i1,a4,i4)') 'broydDa', iDatom, 'Ddir', idir, 'qInd', iqpt

   IF (hybdat%l_calhf) THEN
      OPEN (57,file='hf_broyd',form='unformatted',status='unknown')
   ELSE
      !OPEN (57,file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48),form='unformatted',status='unknown')
      OPEN (57,file=filename,form='unformatted',status='unknown')
   ENDIF

   READ(57) nextIter,alpha,diffDenVec,inDenVec

   CLOSE(57)

END SUBROUTINE readLastIterInAndDiffDen

SUBROUTINE writeLastIterInAndDiffDen(hybdat,vecLen,nextIter,alpha,inDenVec,diffDenVec, idir, iqpt, iDatom)

   TYPE(t_hybdat), INTENT(IN) :: hybdat
   INTEGER,        INTENT(IN) :: nextIter,vecLen
   REAL,           INTENT(IN) :: alpha
   REAL,           INTENT(IN) :: inDenVec(vecLen), diffDenVec(vecLen)
   integer,        intent(in) :: idir
   integer,        intent(in) :: iqpt
   integer,        intent(in) :: iDatom
   character(len=21)           :: filename

   INTEGER             :: mode

   !CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent
   write(filename, '(a7,i1,a4,i1,a4,i4)') 'broydDa', iDatom, 'Ddir', idir, 'qInd', iqpt

   IF (hybdat%l_calhf) THEN
      OPEN (57,file='hf_broyd',form='unformatted',status='unknown')
   ELSE
      !OPEN (57,file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48), form='unformatted',status='unknown')
      OPEN (57,file=filename, form='unformatted',status='unknown')
   ENDIF

   WRITE(57) nextIter,alpha,diffDenVec,inDenVec

   CLOSE(57)

END SUBROUTINE writeLastIterInAndDiffDen

SUBROUTINE readUVec(input,hybdat,vecLen,relIter,currentIter,uVec, idir, iqpt, iDatom)

   TYPE(t_input),  INTENT(IN)  :: input
   TYPE(t_hybdat), INTENT(IN)  :: hybdat
   INTEGER,        INTENT(IN)  :: vecLen, relIter, currentIter
   integer,        intent(in)  :: idir
   integer,        intent(in)  :: iqpt
   integer,        intent(in)  :: iDatom
   REAL,           INTENT(OUT) :: uVec(vecLen)

   character(len=21)           :: filename
   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   !CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   write(filename, '(a7,i1,a4,i1,a4,i4)') 'broydDa', iDatom, 'Ddir', idir, 'qInd', iqpt
   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      !OPEN (59,file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48)//'.'//CHAR(input%imix+48),access='direct',&
      OPEN (59,file=filename//'.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter+relIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   READ(59,rec=2*npos-1) uVec(:vecLen)

   CLOSE(59)

END SUBROUTINE readUVec

SUBROUTINE writeUVec(input,hybdat,vecLen,currentIter,uVec, idir, iqpt, iDatom)

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybdat), INTENT(IN) :: hybdat
   INTEGER,        INTENT(IN) :: vecLen, currentIter
   REAL,           INTENT(IN) :: uVec(vecLen)
   integer,        intent(in) :: idir
   integer,        intent(in) :: iqpt
   integer,        intent(in) :: iDatom
   character(len=21)           :: filename

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   write(filename, '(a7,i1,a4,i1,a4,i4)') 'broydDa', iDatom, 'Ddir', idir, 'qInd', iqpt
   !CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      !OPEN (59,file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48)//'.'//CHAR(input%imix+48),access='direct',&
      OPEN (59,file=filename//'.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   WRITE(59,rec=2*npos-1) uVec(:vecLen)

   CLOSE(59)

END SUBROUTINE writeUVec

SUBROUTINE readVVec(input,hybdat,vecLen,relIter,currentIter,dfivi,vVec, idir, iqpt, iDatom)

   TYPE(t_input),  INTENT(IN)  :: input
   TYPE(t_hybdat), INTENT(IN)  :: hybdat
   INTEGER,        INTENT(IN)  :: vecLen, relIter, currentIter
   integer,        intent(in)  :: idir
   integer,        intent(in)  :: iqpt
   integer,        intent(in)  :: iDatom
   REAL,           INTENT(OUT) :: dfivi
   REAL,           INTENT(OUT) :: vVec(vecLen)

   INTEGER             :: mode, npos
   character(len=21)           :: filename
   INTEGER*8           :: recLen

   !CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   write(filename, '(a7,i1,a4,i1,a4,i4)') 'broydDa', iDatom, 'Ddir', idir, 'qInd', iqpt
   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      !OPEN (59,file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48)//'.'//CHAR(input%imix+48),access='direct',&
      OPEN (59,file=filename//'.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter+relIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   READ(59,rec=2*npos) vVec(:vecLen), dfivi

   CLOSE(59)

END SUBROUTINE readVVec

SUBROUTINE writeVVec(input,hybdat,vecLen,currentIter,dfivi,vVec, idir, iqpt, iDatom)

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybdat), INTENT(IN) :: hybdat
   INTEGER,        INTENT(IN) :: vecLen, currentIter
   REAL,           INTENT(IN) :: dfivi
   REAL,           INTENT(IN) :: vVec(vecLen)
   integer,        intent(in) :: idir
   integer,        intent(in) :: iqpt
   integer,        intent(in) :: iDatom
   character(len=21)           :: filename

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   !CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   write(filename, '(a7,i1,a4,i1,a4,i4)') 'broydDa', iDatom, 'Ddir', idir, 'qInd', iqpt
   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      !OPEN (59,file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48)//'.'//CHAR(input%imix+48),access='direct',&
      OPEN (59,file=filename//'.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   WRITE(59,rec=2*npos) vVec(:vecLen), dfivi

   CLOSE(59)

END SUBROUTINE writeVVec

SUBROUTINE resetBroydenHistory()

   INTEGER           :: i
   LOGICAL           :: l_exist
   CHARACTER(LEN=20) :: filename
   INQUIRE(file='broyd',exist=l_exist)
   IF (l_exist) THEN
      CALL system('rm broyd')
      PRINT *,"Broyden history has been reset."
   END IF
   DO i = 1, 9
      filename = 'broyd.'//CHAR(i+48)
      INQUIRE(file=TRIM(ADJUSTL(filename)),exist=l_exist)
      IF (l_exist) THEN
         CALL system('rm '//TRIM(ADJUSTL(filename)))
      END IF
      filename = 'hf_broyd.'//CHAR(i+48)
      INQUIRE(file=TRIM(ADJUSTL(filename)),exist=l_exist)
      IF (l_exist) THEN
         CALL system('rm '//TRIM(ADJUSTL(filename)))
      END IF
   END DO

END SUBROUTINE resetBroydenHistory

LOGICAL FUNCTION initBroydenHistory(input,hybdat, vecLen, idir, iqpt, iDatom)
! Initializes a Broyden history
! returns true if there already exists a Broyden history

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybdat), INTENT(IN) :: hybdat
   INTEGER,        INTENT(IN) :: vecLen
   integer,        intent(in) :: idir
   integer,        intent(in) :: iqpt
   integer,        intent(in) :: iDatom

   INTEGER*8                  :: recLen
   LOGICAL                    :: l_exist

   character(len=21)           :: filename

   write(filename, '(a7,i1,a4,i1,a4,i4)') 'broydDa', iDatom, 'Ddir', idir, 'qInd', iqpt

   !INQUIRE (file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48)//'.'//CHAR(input%imix+48),exist=l_exist)
   INQUIRE (file=filename//'.'//CHAR(input%imix+48),exist=l_exist)

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      !OPEN (59,file='broydDa'//CHAR(iDatom+48)//'Dd'//CHAR(idir+48)//'q'//CHAR(iqpt+48)//'.'//CHAR(input%imix+48),access='direct',&
      OPEN (59,file=filename//'.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   CLOSE(59)

   initBroydenHistory = l_exist

END FUNCTION initBroydenHistory

END MODULE m_broyd_io
