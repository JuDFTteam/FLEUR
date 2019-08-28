!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_rwinp
      use m_juDFT
      CONTAINS
      SUBROUTINE rw_inp(&
     &                  ch_rw,atoms,obsolete,vacuum,input,stars,sliceplot,banddos,&
     &                  cell,sym,xcpot,noco,oneD,hybrid,kpts,&
     &                  noel,namex,relcor,a1,a2,a3,dtild_opt,name_opt)

!*********************************************************************
!* This subroutine reads or writes an inp - file on unit iofile      *
!* for ch_rw = 'R' read, ch_rw = 'W' write.                          *
!*                                                           Gustav  *
!*********************************************************************
      USE m_calculator
      USE m_types
  
      IMPLICIT NONE
! ..
! ..   Arguments ..
      CHARACTER,INTENT(IN)          :: ch_rw

      TYPE(t_input),INTENT(INOUT)   :: input
      TYPE(t_sym),INTENT(INOUT)     :: sym
      TYPE(t_stars),INTENT(INOUT)   :: stars 
      TYPE(t_atoms),INTENT(INOUT)   :: atoms
      TYPE(t_vacuum),INTENT(INOUT)   :: vacuum
      TYPE(t_obsolete),INTENT(INOUT) :: obsolete
      TYPE(t_kpts),INTENT(INOUT)     :: kpts
      TYPE(t_oneD),INTENT(INOUT)     :: oneD
      TYPE(t_hybrid),INTENT(INOUT)   :: hybrid
      TYPE(t_cell),INTENT(INOUT)     :: cell
      TYPE(t_banddos),INTENT(INOUT)  :: banddos
      TYPE(t_sliceplot),INTENT(INOUT):: sliceplot
      TYPE(t_xcpot_inbuild),INTENT(INOUT)    :: xcpot
      TYPE(t_noco),INTENT(INOUT)     :: noco
    
      REAL,INTENT(INOUT)           :: a1(3),a2(3),a3(3)
      CHARACTER(len=3),INTENT(OUT) :: noel(atoms%ntype)
      CHARACTER(len=4),INTENT(OUT) :: namex 
      CHARACTER(len=12),INTENT(OUT):: relcor
      REAL,INTENT(IN),OPTIONAL     :: dtild_opt
      CHARACTER(len=8),INTENT(IN),OPTIONAL:: name_opt(10)



      CHARACTER(len=8) :: name(10)

!+lda+u
      REAL    u,j
      INTEGER l, i_u
      LOGICAL l_amf
      CHARACTER(len=3) ch_test
      NAMELIST /ldaU/ l,u,j,l_amf
!-lda+u
!+odim
      INTEGER MM,vM,m_cyl
      LOGICAL invs1,zrfs1
      INTEGER chi,rot
      LOGICAL d1,band
      NAMELIST /odim/ d1,MM,vM,m_cyl,chi,rot,invs1,zrfs1
!-odim
! ..
! ..  Local Variables
      REAL     ::scpos  ,zc,dtild   
      INTEGER  ::nw,idsprs
      INTEGER ieq,i,k,na,n,ilo
      REAL s3,ah,a,hs2,rest
      LOGICAL l_hyb,l_sym,ldum
      INTEGER :: ierr, intDummy
! ..
!...  Local Arrays
      CHARACTER :: helpchar(atoms%ntype)
      CHARACTER(len=  4) :: chntype
      CHARACTER(len= 41) :: chform
      CHARACTER(len=100) :: line

!     added for HF and hybrid functionals
      REAL                  ::  aMix,omega
      INTEGER               :: idum
      CHARACTER (len=1)     ::  check

      IF (PRESENT(dtild_opt)) dtild=dtild_opt
      IF (PRESENT(name_opt)) name=name_opt

!     Initialize variables
      l_hyb = .false.

!---------------------------------------------------------------------
      IF (ch_rw.eq.'r') THEN
!--------------------------------------------------------------------
      OPEN (5,file='inp',form='formatted',status='old')
      
      !default not read in in old inp-file
      input%qfix=2
!
      a1(:) = 0
      a2(:) = 0
      a3(:) = 0
      

      WRITE (6,*) '-------- dump of inp-file ------------'
!
      !<-- Added possibility to define variables here

      DO
         READ (UNIT = 5,FMT = 7182,END=77,ERR=77) ch_test
         BACKSPACE(5)
         IF (ch_test   /="def") EXIT
         READ(unit = 5,FMT="(4x,a)") line
         n = INDEX(line,"=")
         IF (n == 0.OR.n>len_TRIM(line)-1) STOP&
     &        "Error in variable definitions"
         CALL ASSIGN_var(line(:n-1),evaluate(line(n+1:)))
      ENDDO

      ! check if isec1 consists of 2 or 3 digits
      READ(UNIT=5,FMT='(29x,a)') check
      BACKSPACE 5
      IF( check .eq. ',' ) THEN
        READ (UNIT=5,FMT=8000,END=99,ERR=99) &
     &                input%strho,input%film,banddos%dos,intDummy,banddos%ndir,input%secvar
        WRITE (6,9000) input%strho,input%film,banddos%dos,99,banddos%ndir,input%secvar
 8000 FORMAT (6x,l1,6x,l1,5x,l1,7x,i2,6x,i2,8x,l1)
      ELSE
        READ (UNIT=5,FMT=8001,END=99,ERR=99) &
     &                input%strho,input%film,banddos%dos,intDummy,banddos%ndir,input%secvar
        WRITE (6,9000) input%strho,input%film,banddos%dos,99,banddos%ndir,input%secvar
 8001 FORMAT (6x,l1,6x,l1,5x,l1,7x,i3,6x,i2,8x,l1)
      END IF

!
      READ (UNIT=5,FMT=7000,END=99,ERR=99) name
      input%comment = name
      WRITE (6,9010) name
 7000 FORMAT (10a8)
!
      READ (UNIT=5,FMT=7020,END=99,ERR=99)&
     &     cell%latnam,sym%namgrp,sym%invs,sym%zrfs,sym%invs2,input%jspins,noco%l_noco
      WRITE (6,9020)&
     &     cell%latnam,sym%namgrp,sym%invs,sym%zrfs,sym%invs2,input%jspins,noco%l_noco
 7020 FORMAT (a3,1x,a4,6x,l1,6x,l1,7x,l1,8x,i1,8x,l1,5x,l1)
!
      IF ((cell%latnam.EQ.'squ').OR.(cell%latnam.EQ.'hex').OR.&
     &    (cell%latnam.EQ.'c-b').OR.(cell%latnam.EQ.'hx3').OR.&
     &    (cell%latnam.EQ.'fcc').OR.(cell%latnam.EQ.'bcc')) THEN
         READ (UNIT = 5,FMT =*,iostat = ierr) a1(1)
         IF (ierr /= 0) THEN
            BACKSPACE(5)
            READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
            a1(1)      = evaluatefirst(line)
         ENDIF
         WRITE (6,9030) a1(1)
      ELSEIF ((cell%latnam.EQ.'c-r').OR.(cell%latnam.EQ.'p-r')) THEN
         READ (UNIT = 5,FMT=*,iostat=ierr) a1(1),a2(2)
         IF (ierr /= 0) THEN
            BACKSPACE(5)
            READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
            a1(1)      = evaluatefirst(line)
            a2(2)      = evaluatefirst(line)
         ENDIF
         WRITE (6,9030) a1(1),a2(2)
      ELSEIF (cell%latnam.EQ.'obl') THEN
         READ (UNIT = 5,FMT =*,iostat= ierr) a1(1),a1(2)
         IF (ierr /= 0) THEN
            BACKSPACE(5)
            READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
            a1(1)      = evaluatefirst(line)
            a1(2)      = evaluatefirst(line)
         ENDIF
         READ (UNIT = 5,FMT =*,iostat= ierr) a2(1),a2(2)
         IF (ierr /= 0) THEN
            BACKSPACE(5)
            READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
            a2(1)      = evaluatefirst(line)
            a2(2)      = evaluatefirst(line)
         ENDIF
         WRITE (6,9030) a1(1),a1(2)
         WRITE (6,9030) a2(1),a2(2)
      ELSEIF (cell%latnam.EQ.'any') THEN
          READ (UNIT=5,FMT=*,iostat= ierr) a1
          IF (ierr /= 0) THEN
             BACKSPACE(5)
             READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
             a1(1)       = evaluatefirst(line)
             a1(2)        = evaluatefirst(line)
             a1(3)        = evaluatefirst(line)
          ENDIF
          READ (UNIT=5,FMT=*,iostat= ierr) a2
          IF (ierr /= 0) THEN
             BACKSPACE(5)
             READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
             a2(1)       = evaluatefirst(line)
             a2(2)        = evaluatefirst(line)
             a2(3)        = evaluatefirst(line)
          ENDIF
          WRITE (6,9030) a1(1),a1(2),a1(3)
          WRITE (6,9030) a2(1),a2(2),a2(3)
      ELSE
          WRITE (6,*) 'rw_inp: cell%latnam ',cell%latnam,' unknown'
           CALL juDFT_error("Unkown lattice name",calledby="rw_inp")
      ENDIF
!
!
      IF (cell%latnam.EQ.'squ') THEN
         a2(2) = a1(1)
      END IF
!
!     Centered rectangular, special case for bcc(110)
!
      IF (cell%latnam.EQ.'c-b') THEN
         a = a1(1)
         hs2 = sqrt(2.)*0.5e0
         a1(1) = a*hs2
         a1(2) = -a*0.5e0
         a2(1) = a*hs2
         a2(2) = a*0.5e0
      END IF
!
!     Centered rectangular, general case
!     on input: a ---> half of long diagonal
!               b ---> half of short diagonal
!
      IF (cell%latnam.EQ.'c-r') THEN
         a1(2) = -a2(2)
         a2(1) =  a1(1)
      END IF
      IF (cell%latnam.EQ.'hex') THEN
         s3 = sqrt(3.)
         ah = a1(1)/2.
         a1(1) = ah*s3
         a1(2) = -ah
         a2(1) = a1(1)
         a2(2) = ah
      END IF
      IF (cell%latnam.EQ.'hx3') THEN
         s3 = sqrt(3.)
         ah = a1(1)/2.
         a1(1) = ah
         a1(2) = -ah*s3
         a2(1) = a1(1)
         a2(2) = -a1(2)
      END IF

      input%scaleA1 = 1.0
      input%scaleA2 = 1.0
      input%scaleC  = 1.0

      IF (sym%namgrp.EQ.'any ') THEN
        INQUIRE (file='sym.out',exist=l_sym)
        IF (.not.l_sym)&
     &       CALL juDFT_error(&
     &       "for namgrp ='any' please provide a sym.out -file !"&
     &       ,calledby ="rw_inp")
      ENDIF
      IF (cell%latnam.EQ.'any') THEN
!        CALL juDFT_error("please specify lattice type (squ,p-r,c-r,hex,hx3,obl)",calledby="rw_inp")
        READ (UNIT=5,FMT=*,iostat=ierr) a3(1),a3(2),a3(3),vacuum%dvac,input%scaleCell
        IF (ierr /= 0) THEN
           BACKSPACE(5)
           READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
           a3(1)           = evaluatefirst(line)
           a3(2)           = evaluatefirst(line)
           a3(3)           = evaluatefirst(line)
           vacuum%dvac     = evaluatefirst(line)
           input%scaleCell = evaluatefirst(line)
        ENDIF
        WRITE (6,9031) a3(1),a3(2),a3(3),vacuum%dvac,input%scaleCell
      ELSE
         READ (UNIT = 5,FMT =*,iostat= ierr) vacuum%dvac,dtild,input%scaleCell
         IF (ierr /= 0) THEN
            BACKSPACE(5)
            READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
            vacuum%dvac     = evaluatefirst(line)
            dtild           = evaluatefirst(line)
            input%scaleCell = evaluatefirst(line)
        ENDIF
        WRITE (6,9030) vacuum%dvac,dtild,input%scaleCell
        a3(3) = dtild
      ENDIF
!
      READ (UNIT=5,FMT=7110,END=99,ERR=99) namex,relcor,aMix,omega
      namex = TRIM(ADJUSTL(namex))
 7110 FORMAT (a4,3x,a12,2f6.3)
      IF ((namex.EQ.'pw91').OR.(namex.EQ.'l91 ').OR.&
     &    (namex.EQ.'pbe') .OR.(namex.EQ.'rpbe').OR.&
     &    (namex.EQ.'Rpbe').OR.(namex.EQ.'wc')  .OR.&
     &    (namex.EQ.'pbe0').OR.(namex.EQ.'hse ').OR.&
     &    (namex.EQ.'lhse').OR.(namex.EQ.'vhse')) THEN                    ! some defaults
         obsolete%lwb=.false. ; obsolete%ndvgrd=6; idsprs=0 ; obsolete%chng=-0.1e-11
      ENDIF
      ! set mixing and screening for variable HSE functional
      WRITE (6,9040) namex,relcor
     
! look what comes in the next two lines
!
      READ (UNIT=5,FMT=7182,END=77,ERR=77) ch_test
      IF (ch_test.EQ.'igr') THEN                          ! GGA input
         BACKSPACE (5)
         READ (UNIT=5,FMT=7121,END=99,ERR=99)&
     &                   idum,obsolete%lwb,obsolete%ndvgrd,idsprs,obsolete%chng
         IF (idsprs.ne.0)&
     &        CALL juDFT_warn("idsprs no longer supported in rw_inp")
         WRITE (6,9121) idum,obsolete%lwb,obsolete%ndvgrd,idsprs,obsolete%chng
 7121    FORMAT (5x,i1,5x,l1,8x,i1,8x,i1,6x,d10.3)

         READ (UNIT=5,FMT=7182,END=77,ERR=77) ch_test
         IF (ch_test.EQ.'igg') THEN                      ! GGA 2nd line
           GOTO 76
         ELSEIF (ch_test.EQ.'&od') THEN                  ! continue with 1D
           GOTO 78
         ELSE
           oneD%odd%d1 = .false.
           GOTO 77
         ENDIF
      ELSEIF (ch_test.EQ.'&od') THEN                  ! continue with 1D
         GOTO 78
      ELSEIF ( ch_test .eq. 'gcu' ) then              ! HF
        BACKSPACE (5)
        READ (UNIT=5,FMT=7999,END=99,ERR=99) hybrid%gcutm1,hybrid%tolerance1,&
     &     hybrid%ewaldlambda,hybrid%lexp,hybrid%bands1
        WRITE (6,9999) hybrid%gcutm1,hybrid%tolerance1,hybrid%ewaldlambda,hybrid%lexp,hybrid%bands1
 7999   FORMAT (6x,f8.5,6x,f10.8,8x,i2,6x,i2,7x,i4)
 9999   FORMAT ('gcutm=',f8.5,',mtol=',f10.8,',lambda=',i2,&
     &          ',lexp=',i2,',bands=',i4)
         oneD%odd%d1 = .false.
         goto 76
      ELSE
         GOTO 77
      ENDIF
!-odim
      READ (UNIT=5,FMT=7182,END=99,ERR=99) ch_test
 7182 FORMAT (a3)
   78 IF (ch_test.EQ.'&od') THEN
        BACKSPACE (5)
        READ (5,odim)
        oneD%odd%d1 = d1 ; oneD%odd%mb = vM ; oneD%odd%M = MM ; oneD%odd%m_cyl = m_cyl
        oneD%odd%chi = chi ; oneD%odd%rot = rot
        oneD%odd%invs = invs1 ; oneD%odd%zrfs = zrfs1
        WRITE (6,8182) d1,MM,vM,m_cyl,chi,rot,invs1,zrfs1
      END IF
!+odim
      GOTO 76
   77 BACKSPACE (5)                                ! continue with atoms
   76 IF (ch_test /= '&od') THEN
        WRITE (6,*) '   '
        oneD%odd%d1 = .false.
        oneD%odd%M = 1 ; oneD%odd%mb = 1 ; oneD%odd%m_cyl = 1
        oneD%odd%chi = 1 ; oneD%odd%rot = 1
        oneD%odd%invs = .FALSE. ; oneD%odd%zrfs = .FALSE.
      END IF
      READ (UNIT=5,FMT=*,END=99,ERR=99) atoms%ntype
      WRITE (6,9050) atoms%ntype
!
      na = 0
      READ (UNIT=5,FMT=7110,END=99,ERR=99)
      WRITE (6,9060)
      atoms%n_u = 0
      atoms%n_hia = 0
      atoms%n_gf = 0
      atoms%n_j0 = 0
      DO n=1,atoms%ntype
!
         READ (UNIT=5,FMT=7140,END=99,ERR=99) noel(n),atoms%nz(n),&
     &                      atoms%ncst(n),atoms%lmax(n),atoms%jri(n),atoms%rmt(n),atoms%dx(n)
         WRITE (6,9070) noel(n),atoms%nz(n),atoms%ncst(n),atoms%lmax(n),atoms%jri(n),&
     &                      atoms%rmt(n),atoms%dx(n)
 7140    FORMAT (a3,i3,3i5,2f10.6)
!
!+lda+u
         READ (UNIT=5,FMT=7180,END=199,ERR=199) ch_test
 7180    FORMAT (a3)
         IF (ch_test.EQ.'&ld') THEN
            l=0 ; u=0.0 ; j=0.0 ; l_amf = .false.
            BACKSPACE (5)
            READ (5,ldaU)
            atoms%n_u = atoms%n_u + 1
            atoms%lda_u(atoms%n_u)%l = l
            atoms%lda_u(atoms%n_u)%u = u
            atoms%lda_u(atoms%n_u)%j = j
            atoms%lda_u(atoms%n_u)%l_amf = l_amf
            atoms%lda_u(atoms%n_u)%atomType = n
            WRITE (6,8180) l,u,j,l_amf
         END IF
 199     CONTINUE
!-lda+u
!
!---> read extra info for local orbitals, and l_geo. if l_geo=T
!---> calculate force on this atom.
!---> p.kurz 97-06-05
!
!     add parameters lcutm and select for HF and hybrid functionals
        IF ( namex=='hf  ' .OR. namex=='pbe0' .OR. namex=='exx '&
     &       .OR. namex=='hse ' .OR. namex=='vhse' ) THEN
          l_hyb = .TRUE.
          READ (UNIT=5,FMT=7160,END=99,ERR=99) atoms%neq(n),&
     &                  atoms%l_geo(n),hybrid%lcutm1(n),hybrid%select1(1,n),hybrid%select1(2,n),&
     &                  hybrid%select1(3,n),hybrid%select1(4,n),atoms%nlo(n),&
     &                  (atoms%llo(ilo,n),ilo=1,atoms%nlo(n))
 7160     FORMAT (i2,8x,l1,7x,i2,8x,i2,1x,i2,1x,i2,1x,i2,5x,i2,5x,60i3)
          WRITE (6,9090) atoms%neq(n),atoms%l_geo(n),hybrid%lcutm1(n),hybrid%select1(1,n),&
     &       hybrid%select1(2,n),hybrid%select1(3,n),hybrid%select1(4,n),atoms%nlo(n),&
     &       (atoms%llo(ilo,n),ilo=1,atoms%nlo(n))
        ELSE
          READ (UNIT=5,FMT=7161,END=99,ERR=99) atoms%neq(n),&
     &                    atoms%l_geo(n),atoms%nlo(n),(atoms%llo(ilo,n),ilo=1,atoms%nlo(n))
 7161     FORMAT (i2,8x,l1,5x,i2,5x,60i3)
          WRITE (6,9091) atoms%neq(n),atoms%l_geo(n),atoms%nlo(n),&
     &                                    (atoms%llo(ilo,n),ilo=1,atoms%nlo(n))
        END IF 
!
         DO ieq=1,atoms%neq(n)
            na = na + 1
            READ (UNIT = 5,FMT = 7170,iostat = ierr)&
     &                      (atoms%taual(i,na),i=1,3),scpos
            IF (ierr.NE.0) THEN
               BACKSPACE(5)
               !<-- read positions with new format
               READ (UNIT = 5,FMT ="(a)",END = 99,ERR = 99) line
               atoms%taual(1,na) = evaluatefirst(line)
               atoms%taual(2,na) = evaluatefirst(line)
               atoms%taual(3,na) = evaluatefirst(line)
               scpos = evaluatefirst(line)
               IF (scpos == 0.0)  scpos          = 1.0
               !>
            ENDIF
            WRITE (6,9100) (atoms%taual(i,na),i=1,3),scpos
 7170       FORMAT (4f10.6)
            IF (scpos.EQ.0.) scpos = 1.
            DO i = 1,2
               atoms%taual(i,na) = atoms%taual(i,na)/scpos
            ENDDO
            IF (.not.input%film) atoms%taual(3,na) = atoms%taual(3,na)/scpos
!+odim
! in 1D case all the coordinates are given cartesian'ly
            IF (oneD%odd%d1) THEN
               atoms%taual(1,na) = atoms%taual(1,na)/a1(1)
               atoms%taual(2,na) = atoms%taual(2,na)/a2(2)
            END IF
!-odim
         ENDDO
         READ (5,*)
         WRITE (6,9060)
      ENDDO
!
      READ (UNIT=5,FMT=7210,END=99,ERR=99) stars%gmax,xcpot%gmaxxc
      WRITE (6,9110) stars%gmax,xcpot%gmaxxc
      stars%gmaxInit = stars%gmax
 7210 FORMAT (2f10.6)
!
      INQUIRE(file='fl7para',exist=ldum)  ! fl7para must not exist for input%gw=2
      IF (input%gw.eq.-1.and.ldum) THEN         ! in the first run of rw_inp
        ldum = .true.                     ! (then, input%gw=-1 at this point).
      ELSE                                !
        ldum = .false.                    !
      ENDIF                               !
      input%gw = 0
      READ (UNIT=5,FMT=7220,END=99,ERR=7215)&
           &                                       input%vchk,input%cdinf,ldum,input%gw,input%gw_neigd
      if (ldum) call judft_error("pot8 not longer supported")
7215  IF(input%strho) input%gw=0
      if (ldum) call judft_error("pot8 not longer supported")
      IF(input%gw.eq.2 .OR. l_hyb) THEN
         IF(ldum)  CALL juDFT_error&
     &        ("Remove fl7para before run with gw = 2!",calledby&
     &        ="rw_inp")
         IF(input%gw_neigd==0)  CALL juDFT_error("No numbands-value given."&
     &        ,calledby ="rw_inp")
      ELSE
        INQUIRE(file='QGpsi',exist=ldum)
        IF(ldum)           CALL juDFT_error&
     &       ("QGpsi exists but gw /= 2 in inp.",calledby ="rw_inp")
      ENDIF
  
      BACKSPACE(5)                                         ! Make sure that input%vchk,input%cdinf,obsolete%pot8 are all given.
      READ (UNIT=5,FMT=7220,END=99,ERR=99) input%vchk,input%cdinf,ldum
      if (ldum) call judft_error("pot8 not longer supported")
      WRITE (6,9120) input%vchk,input%cdinf,.false.,input%gw,input%gw_neigd
 7220 FORMAT (5x,l1,1x,6x,l1,1x,5x,l1,1x,3x,i1,1x,9x,i4)
!
      DO i=1,100 ; line(i:i)=' ' ; ENDDO
      
      input%eig66(2)=.false.
      
      READ (UNIT=5,FMT=6000,END=99,ERR=99)&
     &                idum,ldum,input%l_f,input%eonly,input%eig66(1)!,input%eig66(2)
      WRITE (6,9130) 0,.false.,input%l_f,input%eonly,input%eig66(1)!,input%eig66(2)
 6000 FORMAT (4x,i1,8x,l1,5x,l1,7x,l1,7x,l1)
!
!+roa
      WRITE (chntype,'(i4)') 2*atoms%ntype
      chform = '('//chntype//'i3 )'
      READ (UNIT=5,FMT=chform,END=99,ERR=99) &
     &                (atoms%lnonsph(n),n=1,atoms%ntype)!,(hybrid%lcutwf(n),n=1,atoms%ntype)
      WRITE (6,FMT=chform) (atoms%lnonsph(n),n=1,atoms%ntype)!,(hybrid%lcutwf(n),n=1,atoms%ntype)
 6010 FORMAT (25i3)
!
      READ (UNIT=5,FMT=6010,END=99,ERR=99) nw,obsolete%lepr
      IF (nw.ne.1) CALL juDFT_error("Multiple window calculations not supported")
      WRITE (6,9140) nw,obsolete%lepr
!
      zc=0.0
      READ (UNIT=5,FMT=*,END=99,ERR=99)
         !WRITE (6,'(a8,i2)') 'Window #',nw
!
      READ (UNIT=5,FMT=6040,END=99,ERR=99) input%ellow,input%elup,input%zelec
      WRITE (6,9150) input%ellow,input%elup,input%zelec
6040  FORMAT (4f10.5)
      zc = zc + input%zelec
!
      READ (UNIT=5,FMT='(f10.5)',END=99,ERR=99) input%rkmax
      WRITE (6,FMT='(f10.5,1x,A)') input%rkmax, '=kmax'
     
      READ (UNIT=5,FMT=8010,END=99,ERR=99) input%gauss,input%delgau,input%tria
      WRITE (6,9160) input%gauss,input%delgau,input%tria
 8010 FORMAT (6x,l1,f10.5,5x,l1)
!
      input%tkb=input%delgau
!
      READ(5,fmt='(27x,l1)',END=99,ERR=99) noco%l_soc
      DO i=1,100 ; line(i:i)=' ' ; ENDDO
      BACKSPACE(5)
      READ(5,fmt='(A)',END=99,ERR=99) line
      BACKSPACE(5)
      IF (line(9:10)=='pi') THEN
        READ(5,fmt='(f8.4)') noco%theta
        noco%theta= noco%theta*4.*ATAN(1.)
      ELSE
        READ(5,fmt='(f10.6)',END=99,ERR=99) noco%theta
      ENDIF
      BACKSPACE(5)
      IF (line(19:20)=='pi') THEN
        READ(5,fmt='(10x,f8.4)',END=99,ERR=99) noco%phi
        noco%phi= noco%phi*4.*ATAN(1.)
      ELSE
        READ(5,fmt='(10x,f10.6)',END=99,ERR=99) noco%phi
      ENDIF
      IF ( line(30:34)=='spav=' ) THEN
        BACKSPACE(5)
        READ(5,fmt='(34x,l1)',END=99,ERR=99) noco%l_spav
      ELSE
        noco%l_spav= .false.
      ENDIF
!!$      IF ( line(37:40)=='off=' ) THEN
!!$        BACKSPACE(5)
!!$        chform= '(40x,l1,1x,'//chntype//'a1)'
!!$        CALL judft_error("soc_opt no longer supported")
!!$      ENDIF
    
      READ (UNIT=5,FMT=8050,END=99,ERR=99)&
     &                 input%frcor,sliceplot%slice,input%ctail
      input%coretail_lmax=99
      input%kcrel=0
      BACKSPACE(5)
      READ (UNIT=5,fmt='(A)') line
      input%l_bmt= ( line(52:56)=='bmt=T' ).or.( line(52:56)=='bmt=t' )
      WRITE (6,9170)  input%frcor,sliceplot%slice,input%ctail
 8050 FORMAT (6x,l1,7x,l1,7x,l1,6x,l1,7x,i1,5x,l1,5x,l1)
      
      ! check if itmax consists of 2 or 3 digits
      READ(unit=5,FMT='(8x,a)') check
      BACKSPACE 5

      IF( check .eq. ',' ) THEN
        READ (UNIT=5,FMT=8060,END=99,ERR=99) &
     & input%itmax,input%maxiter,input%imix,input%alpha,input%spinf
        WRITE (6,9180) input%itmax,input%maxiter,input%imix,input%alpha,input%spinf
 8060   FORMAT (6x,i2,9x,i3,6x,i2,7x,f6.2,7x,f6.2)
      ELSE
        READ (UNIT=5,FMT=8061,END=99,ERR=99) &
     & input%itmax,input%maxiter,input%imix,input%alpha,input%spinf
        WRITE (6,9180) input%itmax,input%maxiter,input%imix,input%alpha,input%spinf
 8061   FORMAT (6x,i3,9x,i3,6x,i2,7x,f6.2,7x,f6.2)
      END IF
      
      input%preconditioning_param = 0.0 

      chform = '(5x,l1,'//chntype//'f6.2)'
!      chform = '(5x,l1,23f6.2)'
      READ (UNIT=5,FMT=chform,END=99,ERR=99)&
     &                                   input%swsp, (atoms%bmu(i),i=1,atoms%ntype)
      chform = '(6x,l1,'//chntype//'i3 )'
!      chform = '(6x,l1,23i3 )'
      READ (UNIT=5,FMT=chform,END=99,ERR=99)&
     &                                   input%lflip, (atoms%nflip(i),i=1,atoms%ntype)
!-
      chform = '("swsp=",l1,'//chntype//'f6.2)'
!      chform = '("swsp=",l1,23f6.2)'
      WRITE (6,FMT=chform) input%swsp, (atoms%bmu(i),i=1,atoms%ntype)
      chform = '("lflip=",l1,'//chntype//'i3 )'
!      chform = '("lflip=",l1,23i3 )'
      WRITE (6,FMT=chform) input%lflip, (atoms%nflip(i),i=1,atoms%ntype)
!-roa
!+stm
      READ (UNIT=5,FMT=8075,END=99,ERR=99)&
     &      banddos%vacdos,vacuum%layers,input%integ,vacuum%starcoeff,vacuum%nstars,&
     &      vacuum%locx(1),vacuum%locy(1),vacuum%locx(2),vacuum%locy(2),vacuum%nstm,vacuum%tworkf
      WRITE (6,9210) banddos%vacdos,vacuum%layers,input%integ,vacuum%starcoeff,vacuum%nstars,&
     &      vacuum%locx(1),vacuum%locy(1),vacuum%locx(2),vacuum%locy(2),vacuum%nstm,vacuum%tworkf
 8075 FORMAT (7x,l1,8x,i2,7x,l1,6x,l1,8x,i2,4(4x,f5.2),6x,i1,8x,f10.6)
!-stm
      IF (banddos%vacdos) THEN
        IF (input%integ) THEN
          READ (UNIT=5,FMT=8076,END=99,ERR=99)&
     &                    ((vacuum%izlay(i,k),k=1,2),i=1,vacuum%layers)
          WRITE (6,9220) ((vacuum%izlay(i,k),k=1,2),i=1,vacuum%layers)
 8076     FORMAT (10(2(i3,1x),1x))
        ELSE
          READ (UNIT=5,FMT=8077,END=99,ERR=99)&
     &                    (vacuum%izlay(i,1),i=1,vacuum%layers)
          WRITE (6,9230) (vacuum%izlay(i,1),i=1,vacuum%layers)
 8077     FORMAT (20(i3,1x))
        END IF
      ELSE
        READ (UNIT=5,FMT=*,END=99,ERR=99)
        WRITE (6,fmt='(1x)')
      END IF
!
      band = .false.
      READ (UNIT=5,FMT=8050,END=992,ERR=992) sliceplot%iplot,input%score,sliceplot%plpot,band
      WRITE (6,9240) sliceplot%iplot,input%score,sliceplot%plpot,band
      IF (band) THEN
        banddos%dos=.true. ; banddos%ndir = -4
      ENDIF
      GOTO 993
 992  BACKSPACE(5)
      READ (UNIT=5,FMT=8050,END=99,ERR=99) sliceplot%iplot,input%score,sliceplot%plpot
      WRITE (6,9240) sliceplot%iplot,input%score,sliceplot%plpot,band
!
 993  READ (UNIT=5,FMT='(i3,2f10.6,6x,i3,8x,l1)',END=99,ERR=99)&
     &                sliceplot%kk,sliceplot%e1s,sliceplot%e2s,sliceplot%nnne,input%pallst
      WRITE (6,9250) sliceplot%kk,sliceplot%e1s,sliceplot%e2s,sliceplot%nnne,input%pallst
!
      READ (UNIT=5,FMT=8090,END=99,ERR=99) !
                     !input%xa,input%thetad,input%epsdisp,input%epsforce
      WRITE (6,*) "No relaxation with old input anymore" 
      !input%xa,input%thetad,input%epsdisp,input%epsforce
 8090 FORMAT (3x,f10.5,8x,f10.5,9x,f10.5,10x,f10.5)
!

!+/-odim YM : changed to '70' in the format, sometimes caused probl.
      chform = '(6x,'//chntype//'(3i1,1x))'
      READ (UNIT=5,FMT=chform,END=99,ERR=99)&
     &      ((atoms%relax(i,k),i=1,3),k=1,atoms%ntype)
      chform = '("relax ",'//chntype//'(3i1,1x))'
      WRITE (6,FMT=chform) ((atoms%relax(i,k),i=1,3),k=1,atoms%ntype)

! read dos_params! These will be set automatically if not present!
      banddos%e1_dos=0.0
      banddos%e2_dos=-1.0
      banddos%sig_dos=1e-4
      READ (UNIT=5,FMT='(9x,f10.5,10x,f10.5,9x,f10.5)',&
     &     END=98,ERR=98) banddos%e2_dos,banddos%e1_dos,banddos%sig_dos

      kpts%posScale = 1.0
 
! added for exact-exchange or hybrid functional calculations:
! read in the number of k-points and nx,ny and nz given in the last line
! of the input file,
! we demand that the values given there are consistent with the kpts-file

      IF(namex=='hf  '.OR.namex=='pbe0'.OR.namex=='exx '.OR.namex=='hse '.OR.namex=='vhse'.OR.&
         (banddos%dos.AND.(banddos%ndir == -3))) THEN
         READ (UNIT=5,FMT='(5x,i5,4x,i2,4x,i2,4x,i2)',END=98,ERR=98) idum,kpts%nkpt3(1),kpts%nkpt3(2),kpts%nkpt3(3)

         IF(idum.EQ.0) THEN
            WRITE(*,*) ''
            WRITE(*,*) 'nkpt is set to 0.'
            WRITE(*,*) 'For this fleur mode it has to be larger than 0!'
            WRITE(*,*) ''
            CALL juDFT_error("Invalid declaration of k-point set (1)",calledby="rw_inp")
         END IF
      
         IF( kpts%nkpt3(1)*kpts%nkpt3(2)*kpts%nkpt3(3) .ne. idum ) THEN
            WRITE(*,*) ''
            WRITE(*,*) 'nx*ny*nz is not equal to nkpt.'
            WRITE(*,*) 'For this fleur mode this is required!'
            WRITE(*,*) ''
            CALL juDFT_error("Invalid declaration of k-point set (2)",calledby="rw_inp")
         END IF
      END IF

! for a exx calcuation a second mixed basis set is needed to
! represent the response function, its parameters are read in here

      IF(namex=='exx ') THEN
         CALL judft_error("No EXX calculations in this FLEUR version")
        !READ (UNIT=5,FMT='(7x,f8.5,7x,f10.8,7x,i3)',END=98,ERR=98) hybrid%gcutm2,hybrid%tolerance2,hybrid%bands2

        !DO i=1,atoms%ntype
          !READ (UNIT=5,FMT='(7x,i2,9x,i2,1x,i2,1x,i2,1x,i2)',&
            !END IF=98,ERR=98) hybrid%lcutm2(i),hybrid%select2(1,i),hybrid%select2(2,i),&
            !           hybrid%select2(3,i),hybrid%select2(4,i)
        !END DO
        
        !ALLOCATE( hybrid%l_exxc(maxval(atoms%ncst),atoms%ntype) )
        !DO i=1,atoms%ntype
   !       READ(UNIT=5,FMT='(60(2x,l1))',END=98,ERR=98)(hybrid%l_exxc(k,i),k=1,atoms%ncst(i))
       ! END DO
      END IF

 98   CONTINUE
      WRITE (6,'(a,f10.5,a,f10.5,a,f10.5)')&
     &     'emin_dos=',banddos%e2_dos,',emax_dos=',banddos%e1_dos,',sig_dos=',banddos%sig_dos
      CLOSE (5)




!---------------------------------------------------------------------
      ELSEIF ((ch_rw.eq.'W').OR.(ch_rw.eq.'w'))  THEN
!---------------------------------------------------------------------

      IF (ch_rw.eq.'W') THEN
      OPEN (5,file='inp_new',form='formatted',status='unknown')
      REWIND (5)
      ELSE
      OPEN(5,file='inp',form='formatted',status='unknown')
      ENDIF

      IF (namex.EQ.'hf  ' .OR. namex .EQ. 'exx ' .OR. namex .EQ. 'hse '&
     &  .OR. namex.EQ.'vhse' )&
     &  l_hyb = .true.
      WRITE (5,9000) input%strho,input%film,banddos%dos,99,banddos%ndir,input%secvar
 9000 FORMAT ('strho=',l1,',film=',l1,',dos=',l1,',isec1=',i3,&
     &        ',ndir=',i2,',secvar=',l1)
      WRITE (5,9010) name
 9010 FORMAT (10a8)
      WRITE(5,9020) cell%latnam,sym%namgrp,sym%invs,sym%zrfs,sym%invs2,input%jspins,noco%l_noco
 9020 FORMAT (a3,1x,a4,',invs=',l1,',zrfs=',l1,',invs2=',l1,&
     &       ',jspins=',i1,',l_noco=',l1,',l_J=',l1)
!
      IF (cell%latnam.EQ.'c-b') THEN
         a1(1) = sqrt(2.)* a1(1)
      END IF
      IF (cell%latnam.EQ.'hex') THEN
         s3 = sqrt(3.)
         a1(1) = 2*a1(1)/sqrt(3.)
      END IF
      IF (cell%latnam.EQ.'hx3') THEN
         a1(1) = 2*a1(1)
      END IF
!
      IF ((cell%latnam.EQ.'squ').OR.(cell%latnam.EQ.'hex').OR.&
     &    (cell%latnam.EQ.'c-b').OR.(cell%latnam.EQ.'hx3')) THEN
          WRITE (5,9030) a1(1)
      ELSEIF ((cell%latnam.EQ.'c-r').OR.(cell%latnam.EQ.'p-r')) THEN
          WRITE (5,9030) a1(1),a2(2)
      ELSEIF (cell%latnam.EQ.'obl') THEN
          WRITE (5,9030) a1(1),a1(2)
          WRITE (5,9030) a2(1),a2(2)
      ELSEIF (cell%latnam.EQ.'any') THEN
          WRITE (5,9030) a1(1),a1(2),a1(3)
          WRITE (5,9030) a2(1),a2(2),a2(3)
      ELSE
          WRITE (6,*) 'rw_inp: cell%latnam ',cell%latnam,' unknown'
           CALL juDFT_error("Invalid lattice name",calledby="rw_inp")
      ENDIF
!
      IF (cell%latnam.EQ.'any') THEN
        WRITE (5,9031)  a3(1),a3(2),a3(3),vacuum%dvac,input%scaleCell
        dtild = a3(3)
      ELSE
        WRITE (5,9030) vacuum%dvac,dtild,input%scaleCell
        a3(3) = input%scaleCell * dtild
      ENDIF
 9030 FORMAT (3f15.8)
 9031 FORMAT (5f15.8)
      IF (namex.EQ.'vhse') THEN
        WRITE (5,9041) namex,relcor,0.25,0.11
      ELSE
        WRITE (5,9040) namex,relcor
      ENDIF
 9040 FORMAT (a4,3x,a12)
 9041 FORMAT (a4,3x,a12,2f6.3)
      IF ((namex.EQ.'pw91').OR.(namex.EQ.'l91').OR.&
     &    (namex.eq.'pbe').OR.(namex.eq.'rpbe').OR.&
   &    (namex.EQ.'Rpbe').OR.(namex.eq.'wc') ) THEN
        WRITE (5,FMT=9121) idum,obsolete%lwb,obsolete%ndvgrd,0,obsolete%chng
 9121    FORMAT ('igrd=',i1,',lwb=',l1,',ndvgrd=',i1,',idsprs=',i1,&
     &           ',chng=',d10.3)
      ENDIF
      IF( namex.EQ.'hf  ' .OR. namex .EQ. 'exx ' .OR. namex .EQ. 'hse '&
     &    .OR. namex.EQ.'vhse' ) THEN
        WRITE (5,9999) hybrid%gcutm1,hybrid%tolerance1,hybrid%ewaldlambda,hybrid%lexp,hybrid%bands1
        l_hyb = .true.
      END IF


!-odim
      IF (oneD%odd%d1) THEN
        WRITE (5,8182) oneD%odd%d1,oneD%odd%M,oneD%odd%mb,&
     &                      oneD%odd%m_cyl,oneD%odd%chi,oneD%odd%rot,oneD%odd%invs,oneD%odd%zrfs
 8182   FORMAT ('&odim d1=',l1,',MM=',i3,',vM=',i3,&
     &          ',m_cyl=',i3,',chi=',i3,',rot=',i3,&
     &          ',invs1=',l1,',zrfs1=',l1,'/')
      ELSE
        WRITE (5,*) '   '
      END IF
!+odim
      WRITE (5,9050) atoms%ntype
 9050 FORMAT (i3)
      na = 0
      WRITE (5,9060)
 9060 FORMAT ('**********************************')
      i_u = 1
      DO n=1,atoms%ntype
         WRITE (5,9070) noel(n),atoms%nz(n),atoms%ncst(n),atoms%lmax(n),atoms%jri(n),&
     &                       atoms%rmt(n),atoms%dx(n)
 9070    FORMAT (a3,i3,3i5,2f10.6)
!+lda_u
         IF (i_u.LE.atoms%n_u) THEN
            DO WHILE (atoms%lda_u(i_u)%atomType.LT.n)
               i_u = i_u + 1
               IF (i_u.GE.atoms%n_u) EXIT
            END DO
            IF (atoms%lda_u(i_u)%atomType.EQ.n) THEN
               WRITE (5,8180) atoms%lda_u(i_u)%l,atoms%lda_u(i_u)%u,atoms%lda_u(i_u)%j,atoms%lda_u(i_u)%l_amf
 8180          FORMAT ('&ldaU l=',i1,',u=',f4.2,',j=',f4.2,',l_amf=',l1,'/')
            ELSE
               WRITE (5,*) '   '
            ENDIF
         ELSE
            WRITE (5,*) '   '
         END IF
!-lda_u
        IF ( l_hyb ) THEN
          WRITE (5,9090) atoms%neq(n),atoms%l_geo(n),hybrid%lcutm1(n),hybrid%select1(1,n),&
     &          hybrid%select1(2,n),hybrid%select1(3,n),hybrid%select1(4,n),atoms%nlo(n),&
     &          (atoms%llo(ilo,n),ilo=1,atoms%nlo(n))
 9090     FORMAT ( i2,',force =',l1,',lcutm=',i2,',select=',&
     &        i2,',',i2,';',i2,',',i2,',nlo=',i2,',llo=',60i3 )
        ELSE
          WRITE (5,9091) atoms%neq(n),atoms%l_geo(n),atoms%nlo(n),&
     &          (atoms%llo(ilo,n),ilo=1,atoms%nlo(n))
 9091     FORMAT (i2,',force =',l1,',nlo=',i2,',llo=',60i3)
        END IF
         DO ieq=1,atoms%neq(n)
            na = na + 1

            scpos = 1.0
            DO i = 2,9
               rest = ABS(i*atoms%taual(1,na) - NINT(i*atoms%taual(1,na)) )&
                    + ABS(i*atoms%taual(2,na) - NINT(i*atoms%taual(2,na)) )
               IF (rest.LT.(i*0.000001)) EXIT
            ENDDO
            IF (i.LT.10) scpos = real(i)  ! common factor found (x,y)
            IF (.NOT.input%film) THEN           ! now check z-coordinate
              DO i = 2,9
                rest = ABS(i*atoms%taual(3,na) - NINT(i*atoms%taual(3,na)) )
                IF (rest.LT.(i*scpos*0.000001)) THEN
                  scpos = i*scpos
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            DO i = 1,2
               atoms%taual(i,na) = atoms%taual(i,na)*scpos
            ENDDO
            IF (.NOT.input%film) atoms%taual(3,na) = atoms%taual(3,na)*scpos
            IF (input%film) atoms%taual(3,na) = a3(3)*atoms%taual(3,na)/input%scaleCell
!+odim in 1D case all the coordinates are given in cartesian YM
            IF (oneD%odd%d1) THEN
               atoms%taual(1,na) = atoms%taual(1,na)*a1(1)
               atoms%taual(2,na) = atoms%taual(2,na)*a2(2)
            END IF
!-odim
            WRITE (5,9100) (atoms%taual(i,na),i=1,3),scpos
 9100       FORMAT (4f10.6)
         ENDDO
         WRITE (5,9060)
      ENDDO
      IF ((xcpot%gmaxxc.LE.0).OR.(xcpot%gmaxxc.GT.stars%gmax)) xcpot%gmaxxc=stars%gmax
      WRITE (5,9110) stars%gmax,xcpot%gmaxxc
 9110 FORMAT (2f10.6)
      WRITE (5,9120) input%vchk,input%cdinf,.false.,input%gw,input%gw_neigd
 9120 FORMAT ('vchk=',l1,',cdinf=',l1,',pot8=',l1,',gw=',i1,&
     &        ',numbands=',i4)
      WRITE (5,9130) 0,.false.,input%l_f,input%eonly
 9130 FORMAT ('lpr=',i1,',form66=',l1,',l_f=',l1,',eonly=',l1)
      IF ( l_hyb ) THEN
        WRITE (chntype,'(i3)') 2*atoms%ntype
        chform = '('//chntype//'i3 )'
        WRITE (5,FMT=chform) &
     &        (atoms%lnonsph(n),n=1,atoms%ntype),(hybrid%lcutwf(n),n=1,atoms%ntype)
      ELSE
         WRITE (chntype,'(i3)') atoms%ntype
         chform = '('//chntype//'i3 )'
         WRITE (5,FMT=chform) (atoms%lnonsph(n),n=1,atoms%ntype)
      END IF
 9140 FORMAT (25i3)
      WRITE (5,9140) 1,obsolete%lepr

      WRITE (5,'(a)') 'ellow, elup, valence electrons:'
      
      WRITE (5,9150) input%ellow,input%elup,input%zelec
9150  FORMAT (4f10.5)
      WRITE (5,fmt='(f10.5,1x,A)') input%rkmax, '=kmax'
      WRITE (5,9160) input%gauss,input%delgau,input%tria
 9160 FORMAT ('gauss=',l1,f10.5,'tria=',l1)
      WRITE (5,9170) input%frcor,sliceplot%slice,input%ctail
 9170 FORMAT ('frcor=',l1,',slice=',l1,',ctail=',l1)
      WRITE (5,9180) input%itmax,input%maxiter,input%imix,input%alpha,input%spinf
 9180 FORMAT ('itmax=',i3,',maxiter=',i3,',imix=',i2,',alpha=',&
     &        f6.2,',spinf=',f6.2)
!+roa
      WRITE (chntype,'(i3)') atoms%ntype
      chform = '("swsp=",l1,'//chntype//'f6.2)'
      WRITE (5,FMT=chform) input%swsp, (atoms%bmu(i),i=1,atoms%ntype)
      chform = '("lflip=",l1,'//chntype//'i3 )'
      WRITE (5,FMT=chform) input%lflip, (atoms%nflip(i),i=1,atoms%ntype)
!-roa
!+stm
      WRITE (5,9210) banddos%vacdos,vacuum%layers,input%integ,vacuum%starcoeff,vacuum%nstars,&
     &      vacuum%locx(1),vacuum%locy(1),vacuum%locx(2),vacuum%locy(2),vacuum%nstm,vacuum%tworkf
 9210 FORMAT ('vacdos=',l1,',layers=',i2,',integ=',l1,',star=',l1,&
     & ',nstars=',i2,4(4x,f5.2),',nstm=',i1,',tworkf=',f10.6)
!-stm
      IF (banddos%vacdos) THEN
        IF (input%integ) THEN
          WRITE (5,9220) ((vacuum%izlay(i,k),k=1,2),i=1,vacuum%layers)
 9220     FORMAT (10(2(i3,1x),1x))
        ELSE
          WRITE (5,9230) (vacuum%izlay(i,1),i=1,vacuum%layers)
 9230     FORMAT (20(i3,1x))
        END IF
      ELSE
        WRITE (5,*)
      END IF
      band = .false.
      WRITE (5,9240) sliceplot%iplot,input%score,sliceplot%plpot,band
 9240 FORMAT ('iplot=',l1,',score=',l1,',plpot=',l1,',band=',l1)
      WRITE (5,9250) sliceplot%kk,sliceplot%e1s,sliceplot%e2s,sliceplot%nnne,input%pallst
 9250 FORMAT (i3,2f10.6,',nnne=',i3,',pallst=',l1)
      WRITE(5,*) "No relaxation with old input anymore"
      !WRITE (5,9260) input%xa,input%thetad,input%epsdisp,input%epsforce
 9260 FORMAT ('xa=',f10.5,',thetad=',f10.5,',epsdisp=',f10.5,&
     &        ',epsforce=',f10.5)
!+/-gb
      chform = '("relax ",'//chntype//'(3i1,1x))'
      WRITE (5,FMT=chform) ((atoms%relax(i,k),i=1,3),k=1,atoms%ntype)
      WRITE (5,'(a,f10.5,a,f10.5,a,f10.5)')&
     &     'emin_dos=',banddos%e2_dos,',emax_dos=',banddos%e1_dos,',sig_dos=',banddos%sig_dos

      IF (ch_rw.eq.'W') CLOSE (5)
      ELSE
        WRITE (6,*) 'specify either W to write or R to read!'
      ENDIF

      RETURN

      !Error handling only here
  99  WRITE (6,*) 'Error reading inp-file'
      CLOSE (6)
      CALL juDFT_error("error reading inp-file",calledby="rw_inp")

      END SUBROUTINE rw_inp
      END MODULE m_rwinp
