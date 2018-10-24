      MODULE m_firstglance
      use m_juDFT
c
c reads the part of the input file that is necessary to call rw_inp
c
      CONTAINS
      SUBROUTINE first_glance(
     <                        ntype,nop,nat,nlod,layerd,itmax,
     <                        l_kpts,l_qpts,l_gamma,nkpt,nmop,
     <                        nmopq)

      USE m_symdata , ONLY : nammap,ord2,l_c2
      IMPLICIT NONE

      INTEGER,INTENT (OUT) :: ntype,nop,nat,nlod,layerd,itmax
      INTEGER,INTENT (OUT) :: nkpt,nmop(3),nmopq(3)
      LOGICAL,INTENT (OUT) :: l_kpts,l_qpts,l_gamma

!mod: INTEGER ord2(25)               ! Number of 2D symmetry operations
!mod: LOGICAL l_c2(25)               ! whether plane group contains the c_2
!mod: CHARACTER(len=4) :: nammap(20) ! names as in the inp-file

      INTEGER neq,n,na,nlo,line,i,n2spg,nqpt
      CHARACTER*4 namgrp,namex
      CHARACTER*3 latnam
      LOGICAL invs,zrfs

      l_kpts = .false.
      l_qpts=.true.

      OPEN (5,file='inp',form='formatted',status='old')
      
      !<-- skip first lines with definitions
      DO
         READ(5,*) latnam
         if (latnam/="def") exit
      ENDDO
      backspace(5)
      !>

      line = 0
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,8010) latnam,namgrp,invs,zrfs
      line = line + 1
 8010 FORMAT (a3,1x,a4,6x,l1,6x,l1)

      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      IF ((latnam.EQ.'any').OR.(latnam.EQ.'obl')) THEN
        READ (5,*)
        line = line + 1
      ENDIF
      READ (5,'(a4)') namex
      line = line + 1

      READ (5,*,END=77,ERR=77) ntype
      GOTO 79
   77 READ (5,*,END=78,ERR=78) ntype
      GOTO 79
   78 READ (5,*,END=99,ERR=99) ntype
   79 line = line + 1
      nlod = 0
      nat = 0
      DO n = 1,ntype
         READ (5,*)
         line = line + 1
         READ (5,*)
         line = line + 1
         READ (5,*)
         line = line + 1
         IF ( namex=='hf  ' .OR. namex=='pbe0' .OR. namex=='exx '
     +        .OR. namex=='hse ' .OR. namex=='vhse' ) THEN
          READ (5,'(i2,42x,i2)',END=99,ERR=99) neq,nlo
         ELSE
          READ (5,'(i2,14x,i2)',END=99,ERR=99) neq,nlo
         ENDIF
         line = line + 1
         nlod = max(nlo,nlod)
         DO na = 1,neq
            READ (5,*)
            line = line + 1
         ENDDO
         nat = nat + neq
      ENDDO
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      DO n = 1,1!nwdd
         READ (5,*)
         line = line + 1
         READ (5,*)
         line = line + 1
         READ (5,*)
         line = line + 1
      ENDDO
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,'(6x,i2)',END=99,ERR=99) itmax
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,*)
      line = line + 1
      READ (5,'(16x,i2)',END=99,ERR=99) layerd
      line = line + 1

c if (.not.l_kpts) then the kpts-file is missing and we have
c to find out how many k-points to generate...

      INQUIRE (file='QGpsi',exist=l_kpts)
      IF (.not.l_kpts) INQUIRE (file='kpts',exist=l_kpts)
      IF (.not.l_kpts) THEN
        WRITE (6,*) 'No kpts-file exists, trying to generate it'
        DO line = 1,6
          READ (5,*,END=95,ERR=95)
        ENDDO
        READ (5,'(5x,i5,3(4x,i2),7x,l1)',END=97,ERR=97) 
     +                   nkpt,nmop(1),nmop(2),nmop(3),l_gamma
        GOTO 96
 97     BACKSPACE (5)
        READ (5,'(5x,i5,3(4x,i2))',END=98,ERR=98)
     +                       nkpt,nmop(1),nmop(2),nmop(3)
        l_gamma=.false.
        GOTO 96
 98     BACKSPACE (5)
        READ (5,'(5x,i5)',END=95,ERR=95) nkpt
        nmop(1) = 0 ; nmop(2) = 0 ; nmop(3) = 0 ; l_gamma=.false.
        GOTO 96

 95     WRITE (6,*) 'Since you did not provide a kpts-file, you'
        WRITE (6,*) 'should give k-mesh information at the end of'
        WRITE (6,*) 'the inp-file, at least how many k-points, e.g.'
        WRITE (6,*) 'nkpt=  100  -- or give the divisions n in xyz:'
        WRITE (6,*) 'nkpt=   36,nx= 6,ny= 6,nz= 8,gamma=F'
        CLOSE (6)
        STOP
      ENDIF
 96   CONTINUE


! determine the number of symmetry operations
      IF (namgrp.EQ.'any ') THEN
         nop = 48
         CLOSE (5)
         RETURN
      ENDIF
      n2spg = 0
      DO i = 1, 20
        IF (namgrp.EQ.nammap(i)) n2spg = i
      ENDDO
      IF (n2spg == 0 ) THEN
        WRITE (*,*) 'Spacegroup ',namgrp,' not known! Choose one of:'
        WRITE (*,'(20(a4,1x))') (nammap(i),i=1,20)
        CALL juDFT_error("Could not determine spacegroup!",calledby
     +       ="first_glance")
      ENDIF
      IF ( (n2spg.GE.13).AND.(n2spg.LE.17) ) THEN
        IF ( .not.((latnam.EQ.'hx3').OR.(latnam.EQ.'hex')) ) THEN
           CALL juDFT_error
     +          ("Use only hex or hx3 with p3, p3m1, p31m, p6 or p6m!"
     +          ,calledby ="first_glance")
        ENDIF
      ENDIF

      nop = ord2(n2spg)
      IF (invs) THEN
         nop = 2*nop
         IF ( zrfs.and.(.not.l_c2(n2spg)) ) nop = 2*nop
      ELSE
         IF (zrfs) nop = 2*nop
      ENDIF

      CLOSE (5)
      RETURN 
!
! Error
!
  99  WRITE (6,*) 'Error glancing at inp file at line',line
      CLOSE (6)
      CALL juDFT_error("Error glancing at inp file",calledby
     +       ="first_glance")

      END SUBROUTINE first_glance
      END MODULE m_firstglance
