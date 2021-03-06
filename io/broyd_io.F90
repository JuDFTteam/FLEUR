!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_broyd_io

USE m_types
USE m_cdn_io

IMPLICIT NONE

CONTAINS

SUBROUTINE readLastIterInAndDiffDen(hybinp,vecLen,nextIter,alpha,inDenVec,diffDenVec,inDenVecMet,diffDenVecMet)

   TYPE(t_hybinp), INTENT(IN)  :: hybinp
   INTEGER,        INTENT(IN)  :: vecLen
   INTEGER,        INTENT(OUT) :: nextIter
   REAL,           INTENT(OUT) :: alpha
   REAL,           INTENT(OUT) :: inDenVec(vecLen), diffDenVec(vecLen)
   REAL,OPTIONAL,  INTENT(OUT) :: inDenVecMet(vecLen), diffDenVecMet(vecLen)

   INTEGER             :: mode

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   IF (hybdat%l_calhf) THEN
      OPEN (57,file='hf_broyd',form='unformatted',status='unknown')
   ELSE
      OPEN (57,file='broyd',form='unformatted',status='unknown')
   ENDIF

   IF(PRESENT(inDenVecMet)) THEN
      READ(57) nextIter,alpha,diffDenVec,inDenVec,inDenVecMet,diffDenVecMet
   ELSE
      READ(57) nextIter,alpha,diffDenVec,inDenVec
   END IF

   CLOSE(57)

END SUBROUTINE readLastIterInAndDiffDen

SUBROUTINE writeLastIterInAndDiffDen(hybinp,vecLen,nextIter,alpha,inDenVec,diffDenVec,inDenVecMet,diffDenVecMet)

   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: nextIter,vecLen
   REAL,           INTENT(IN) :: alpha
   REAL,           INTENT(IN) :: inDenVec(vecLen), diffDenVec(vecLen)
   REAL,OPTIONAL,  INTENT(IN) :: inDenVecMet(vecLen), diffDenVecMet(vecLen)

   INTEGER             :: mode

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   IF (hybdat%l_calhf) THEN
      OPEN (57,file='hf_broyd',form='unformatted',status='unknown')
   ELSE
      OPEN (57,file='broyd',form='unformatted',status='unknown')
   ENDIF

   IF(PRESENT(inDenVecMet)) THEN
      WRITE(57) nextIter,alpha,diffDenVec,inDenVec,inDenVecMet,diffDenVecMet
   ELSE
      WRITE(57) nextIter,alpha,diffDenVec,inDenVec
   END IF

   CLOSE(57)

END SUBROUTINE writeLastIterInAndDiffDen

SUBROUTINE readUVec(input,hybinp,vecLen,relIter,currentIter,uVec)

   TYPE(t_input),  INTENT(IN)  :: input
   TYPE(t_hybinp), INTENT(IN)  :: hybinp
   INTEGER,        INTENT(IN)  :: vecLen, relIter, currentIter
   REAL,           INTENT(OUT) :: uVec(vecLen)

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter+relIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   READ(59,rec=2*npos-1) uVec(:vecLen)

   CLOSE(59)

END SUBROUTINE readUVec

SUBROUTINE writeUVec(input,hybinp,vecLen,currentIter,uVec)

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: vecLen, currentIter
   REAL,           INTENT(IN) :: uVec(vecLen)

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   WRITE(59,rec=2*npos-1) uVec(:vecLen)

   CLOSE(59)

END SUBROUTINE writeUVec

SUBROUTINE readVVec(input,hybinp,vecLen,relIter,currentIter,dfivi,vVec)

   TYPE(t_input),  INTENT(IN)  :: input
   TYPE(t_hybinp), INTENT(IN)  :: hybinp
   INTEGER,        INTENT(IN)  :: vecLen, relIter, currentIter
   REAL,           INTENT(OUT) :: dfivi
   REAL,           INTENT(OUT) :: vVec(vecLen)

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter+relIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   READ(59,rec=2*npos) vVec(:vecLen), dfivi

   CLOSE(59)

END SUBROUTINE readVVec

SUBROUTINE writeVVec(input,hybinp,vecLen,currentIter,dfivi,vVec)

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: vecLen, currentIter
   REAL,           INTENT(IN) :: dfivi
   REAL,           INTENT(IN) :: vVec(vecLen)

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter-1
   IF (currentIter.GT.input%maxiter+1) npos = MOD(currentIter-2,input%maxiter)+1

   WRITE(59,rec=2*npos) vVec(:vecLen), dfivi

   CLOSE(59)

END SUBROUTINE writeVVec









SUBROUTINE readDeltaNVec(input,hybinp,vecLen,relIter,currentIter,deltaNVec)

   TYPE(t_input),  INTENT(IN)  :: input
   TYPE(t_hybinp), INTENT(IN)  :: hybinp
   INTEGER,        INTENT(IN)  :: vecLen, relIter, currentIter
   REAL,           INTENT(OUT) :: deltaNVec(vecLen)

   INTEGER             :: mode, npos, iter
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd_DN',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd_DN',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   iter = currentIter+relIter
   npos=iter-1
   IF (iter.GT.1) npos = MOD(iter-2,input%maxiter)+1

!   WRITE(1400,*) 'readDeltaNVec: iter,npos: ', iter, npos

   READ(59,rec=npos) deltaNVec(:vecLen)

   CLOSE(59)

END SUBROUTINE readDeltaNVec

SUBROUTINE writeDeltaNVec(input,hybinp,vecLen,currentIter,deltaNVec)

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: vecLen, currentIter
   REAL,           INTENT(IN) :: deltaNVec(vecLen)

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd_DN',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd_DN',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter-1
   IF (currentIter.GT.1) npos = MOD(currentIter-2,input%maxiter)+1

!   WRITE(1500,*) 'writeDeltaNVec: iter,npos: ', currentIter, npos

   WRITE(59,rec=npos) deltaNVec(:vecLen)

   CLOSE(59)

END SUBROUTINE writeDeltaNVec

SUBROUTINE readDeltaFVec(input,hybinp,vecLen,relIter,currentIter,deltaFVec)

   TYPE(t_input),  INTENT(IN)  :: input
   TYPE(t_hybinp), INTENT(IN)  :: hybinp
   INTEGER,        INTENT(IN)  :: vecLen, relIter, currentIter
   REAL,           INTENT(OUT) :: deltaFVec(vecLen)

   INTEGER             :: mode, npos, iter
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd_DF',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd_DF',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   iter = currentIter+relIter
   npos=iter-1
   IF (iter.GT.1) npos = MOD(iter-2,input%maxiter)+1

!   WRITE(1400,*) 'readDeltaFVec: iter,npos: ', iter, npos

   READ(59,rec=npos) deltaFVec(:vecLen)

   CLOSE(59)

END SUBROUTINE readDeltaFVec

SUBROUTINE writeDeltaFVec(input,hybinp,vecLen,currentIter,deltaFVec)

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: vecLen, currentIter
   REAL,           INTENT(IN) :: deltaFVec(vecLen)

   INTEGER             :: mode, npos
   INTEGER*8           :: recLen

   CALL getIOMode(mode)

   ! At the moment broyden IO is mode independent

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd_DF',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd_DF',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter-1
   IF (currentIter.GT.1) npos = MOD(currentIter-2,input%maxiter)+1

!   WRITE(1500,*) 'writeDeltaFVec: iter,npos: ', currentIter, npos

   WRITE(59,rec=npos) deltaFVec(:vecLen)

   CLOSE(59)

END SUBROUTINE writeDeltaFVec


SUBROUTINE writeBroydenOverlapExt(input,hybinp,currentIter,historyLength,&
                                  dNdNLast,dFdFLast,dNdFLast,dFdNLast)

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: currentIter, historyLength
   REAL,           INTENT(IN) :: dNdNLast(input%maxIter)
   REAL,           INTENT(IN) :: dFdFLast(input%maxIter)
   REAL,           INTENT(IN) :: dNdFLast(input%maxIter)
   REAL,           INTENT(IN) :: dFdNLast(input%maxIter)

   INTEGER                    :: recLen, npos

   recLen = 8*4*input%maxIter    ! sizeOfReal*numberOfVectors*vectorLength
   recLen = recLen + 2*8         ! storage for currentIter, historyLength

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broydOvlp',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broydOvlp',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   npos=currentIter-1
   IF (currentIter.GT.1) npos = MOD(currentIter-2,input%maxiter)+1

!   WRITE(1500,*) 'writeBroydenOverlapExt: iter, npos: ', currentIter, npos

   WRITE(59,rec=npos) currentIter, historyLength,&
                      dNdNLast(:input%maxIter), dFdFLast(:input%maxIter), &
                      dNdFLast(:input%maxIter), dFdNLast(:input%maxIter)

   CLOSE(59)

END SUBROUTINE writeBroydenOverlapExt

SUBROUTINE readBroydenOverlaps(input,hybinp,currentIter,historyLength,&
                               dNdNMat,dFdFMat,dNdFMat,dFdNMat)

   TYPE(t_input),  INTENT(IN)    :: input
   TYPE(t_hybinp), INTENT(IN)    :: hybinp
   INTEGER,        INTENT(IN)    :: currentIter, historyLength
   REAL,           INTENT(INOUT) :: dNdNMat(historyLength,historyLength)
   REAL,           INTENT(INOUT) :: dFdFMat(historyLength,historyLength)
   REAL,           INTENT(INOUT) :: dNdFMat(historyLength,historyLength)
   REAL,           INTENT(INOUT) :: dFdNMat(historyLength,historyLength)

   INTEGER                    :: i, j
   INTEGER                    :: recLen, npos, iter
   INTEGER                    :: recIter, recHistLen

   REAL                       :: dNdNLast(input%maxIter)
   REAL                       :: dFdFLast(input%maxIter)
   REAL                       :: dNdFLast(input%maxIter)
   REAL                       :: dFdNLast(input%maxIter)

   recLen = 8*4*input%maxIter    ! sizeOfReal*numberOfVectors*vectorLength
   recLen = recLen + 2*8         ! storage for currentIter, historyLength

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broydOvlp',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broydOvlp',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   dNdNMat = 0.0
   dFdFMat = 0.0
   dNdFMat = 0.0
   dFdNMat = 0.0

   DO i = 1, historyLength

      iter = currentIter - historyLength + i
      npos=iter-1
      IF (iter.GT.1) npos = MOD(iter-2,input%maxiter)+1

!      WRITE(1400,*) 'readBroydenOverlaps: iter,npos: ', iter, npos

      READ(59,rec=npos) recIter, recHistLen,&
                        dNdNLast(:input%maxIter), dFdFLast(:input%maxIter), &
                        dNdFLast(:input%maxIter), dFdNLast(:input%maxIter)

      DO j = 1, recHistLen
         IF ((j-recHistLen+i).LE.0) CYCLE
         dNdNMat(j-recHistLen+i,i) = dNdNLast(j)
         dFdFMat(j-recHistLen+i,i) = dFdFLast(j)
         dNdFMat(j-recHistLen+i,i) = dNdFLast(j)
         dFdNMat(j-recHistLen+i,i) = dFdNLast(j)
      END DO
   END DO

   CLOSE(59)

END SUBROUTINE readBroydenOverlaps







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
   INQUIRE(file='broydOvlp',exist=l_exist)
   IF (l_exist) CALL system('rm broydOvlp')
   INQUIRE(file='hf_broydOvlp',exist=l_exist)
   IF (l_exist) CALL system('rm hf_broydOvlp')

   INQUIRE(file='broyd_DF',exist=l_exist)
   IF (l_exist) CALL system('rm broyd_DF')
   INQUIRE(file='hf_broyd_DF',exist=l_exist)
   IF (l_exist) CALL system('rm hf_broyd_DF')

   INQUIRE(file='broyd_DN',exist=l_exist)
   IF (l_exist) CALL system('rm broyd_DN')
   INQUIRE(file='hf_broyd_DN',exist=l_exist)
   IF (l_exist) CALL system('rm hf_broyd_DN')

END SUBROUTINE resetBroydenHistory

LOGICAL FUNCTION initBroydenHistory(input,hybinp, vecLen)
! Initializes a Broyden history
! returns true if there already exists a Broyden history

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: vecLen

   INTEGER*8                  :: recLen
   LOGICAL                    :: l_exist

   INQUIRE (file='broyd.'//CHAR(input%imix+48),exist=l_exist)

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd.'//CHAR(input%imix+48),access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   CLOSE(59)

   initBroydenHistory = l_exist

END FUNCTION initBroydenHistory

LOGICAL FUNCTION initBroydenHistory2(input,hybinp, vecLen)
! Initializes a Broyden history
! returns true if there already exists a Broyden history

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_hybinp), INTENT(IN) :: hybinp
   INTEGER,        INTENT(IN) :: vecLen

   INTEGER*8                  :: recLen
   LOGICAL                    :: l_exist

   INQUIRE (file='broyd_DF',exist=l_exist)

   recLen=(vecLen+1)*8

   IF (hybdat%l_calhf) THEN
      OPEN (59,file='hf_broyd_DF',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
      OPEN (60,file='hf_broyd_DN',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ELSE
      OPEN (59,file='broyd_DF',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
      OPEN (60,file='broyd_DN',access='direct',&
            recl=recLen,form='unformatted',status='unknown')
   ENDIF

   CLOSE(59)
   CLOSE(60)

   initBroydenHistory2 = l_exist

END FUNCTION initBroydenHistory2

END MODULE m_broyd_io
