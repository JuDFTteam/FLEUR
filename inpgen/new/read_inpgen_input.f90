!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_read_inpgen_input
  USE m_judft
  IMPLICIT NONE
CONTAINS

  SUBROUTINE read_inpgen_input(atom_pos,atom_id,atom_label,amat,div,namex,relcor,dtild&
       input,sym,noco,vacuum,stars,kpts,xcpot,&
       filename)
    !Subroutine reads the old-style input for inpgen
    !if no filename is given std-in is used
    
    REAL,    ALLOCATABLE,INTENT(OUT) :: atom_pos(:, :),atom_id(:)
    CHARACTER(len=20), ALLOCATABLE,INTENT(OUT) :: atom_Label(:)
    REAL,INTENT(out)    :: amat(3,3)
    INTEGER,INTENT(out)::div(3)
    CHARACTER(len=4),INTENT(out)  :: namex
    CHARACTER(len=12),INTENT(out) :: relcor
    REAL,INTENT(OUT) :: dtild

    TYPE(t_input),INTENT(out)   :: input
    TYPE(t_sym),INTENT(OUT)     :: sym
    TYPE(t_noco),INTENT(OUT)    :: noco
    TYPE(t_vacuum),INTENT(OUT)  :: vacuum
    TYPE(t_stars),INTENT(OUT)   :: stars
    TYPE(t_ktps),INTENT(OUT)    :: kpts
    TYPE(t_xcpot),INTENT(OUT)   :: xcpot

    CHARACTER(len=*),INTENT(IN),OPTIONAL::filename
    
    !local variables
    TYPE(t_obsolete)     :: obsolete
    CHARACTER(len=16384) :: buffer
    INTEGER,PARAMTER     :: natmax=9999
    CHARACTER(len=80)    :: title
    LOGICAL              :: cal_symm,checkinp,inistop,cartesian,l_hyb,oldfleur
    INTEGER              :: ngen,i_c,natin,nline,i,j
    INTEGER, ALLOCATABLE :: mmrot(:,:,:)
    REAL,    ALLOCATABLE :: ttr(:, :),atompos(:, :),atomid(:)
    CHARACTER(len=20), ALLOCATABLE :: atomLabel(:)
    REAL :: a1(3),a2(3),a3(3),SCALE(3),factor(3),aa  !lattice definition
    INTEGER:: infh,errfh,warnfh,dbgfh,outfh,symfh    !file handles

    filename=juDFT_string_for_argument("-f")
    INQUIRE(file=filename,exist=l_exist)
    IF (.NOT.l_exist) CALL judft_error("Input file specified is not readable")

    OPEN(97,file=filename)
    OPEN(98,status='scratch')

    CALL  normalize_file(97,98)

    REWIND(98)
    READ(98,"(a)",iostat=ios) input%comment
    namelist_ok=.TRUE.
    DO WHILE(ios.NE.0)
       READ(98,"(a)",iostat=ios) line
       IF (ios.NE.0) EXIT
       IF (line(1:1)=="&") THEN
          !process the namelist
          SELECT CASE(line(2:5)) !e.g. atom
          CASE ('latt')
             CALL process_lattice(line,cell)???
          CASE('inpu')
             CALL process_input(line,input%film,sym%symor)
          CASE('qss ')
             CALL process_qss(line,noco)
          CASE('soc ')
             CALL process_soc(line,noco)
          CASE('shif')
             CALL process_shift(line,atom_pos)
          CASE('fact')
             CALL process_factor(line,atom_pos)
          CASE('exco')
             CALL process_exco(line,xcpot)
          CASE('comp')
             CALL process_comp(line,input%jspins,input%frcor,input%ctail,input%kcrel,stars%gmax,xcpot%gmaxxc,input%rkmax)
          CASE('kpt ')
             CALL process_kpt(line,???)???
          CASE('film')
             CALL process_film(line,vacuum%dvac,vacuum%dtild)???
          CASE('gen ','sym ')
             CALL judft_error("Specifying the symmetries no longer supported in inpgen")
          CASE default
             CALL judft_error(("Unkown input in:"//line))
          END SELECT
       ELSE
          IF (SUM(ABS(cell%amat))>0) THEN
             !cell was set already, so list of atoms follow
             READ(line,*,iostat=ios) n
             IF (ios.NE.0) CALL judft_error(("Surprising error in reading input:"//line))
             ALLOCATE(atom_pos(3,n),atom_label(n),atom_id(n))
             DO i=1,n
                READ(98,"(a)",iostat=ios) line
                IF (ios.NE.0) CALL judft_error(("List of atoms not complete:"//line))
                atom_id(i)=evaluatefirst(line)
                atom_pos(1,i)=evaluatefirst(line)
                atom_pos(2,i)=evaluatefirst(line)
                atom_pos(3,i)=evaluatefirst(line)
                IF(TRIM(ADJUSTL(line)).NE.'') THEN
                   atom_Label(i) = TRIM(ADJUSTL(line))
                ELSE
                   WRITE(atom_Label(i),'(i0)') n
                END IF
             END DO
          ELSE
             !the bravais matrix has to follow
             ???
          ENDIF
       ENDIF
    END DO

    IF (.NOT.ALLOCATED(atompos).OR.SUM(ABS(cell%amat))==0.0) CALL judft_error("input not complete")
    
  END SUBROUTINE read_inpgen_input

    SUBROUTINE process_input(line,film,symor,hybrid)
      CHARACTER(len=*),INTENT(in)::line
      LOGICAL,INTENT(out)::film,symor,hybrid

      INTEGER :: ios
      LOGICAL :: cartesian, cal_symm, checkinp,inistop,oldfleur
      cartesian=.FALSE.
      cal_sym=.FALSE.
      oldfleur=.FALSE.
      NAMELIST /input/ film, cartesian, cal_symm, checkinp, inistop,
     &                 symor, oldfleur, hybrid
     READ(line,input,iostat=ios)
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
     IF (ANY([cal_symm, checkinp,oldfleur])) CALL judft_error("Switches cal_symm, checkinp,oldfleur no longer supported")
   END SUBROUTINE process_input
   SUBROUTINE process_qss(line,noco)
     CHARACTER(len=*),INTENT(in)::line
     TYPE(t_noco),INTENT(INOUT) :: noco
     CHARACTER(len=1000) :: buf

     buf=ADJUSTL(line(5:len_TRIM(line)-1)
          
     READ(line,*,iostat=ios) noco%qss
     noco%l_ss=.TRUE.
     noco%l_noco=.TRUE.
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
   END SUBROUTINE process_qss
   SUBROUTINE process_soc(line,noco)
     CHARACTER(len=*),INTENT(in)::line
     TYPE(t_noco),INTENT(INOUT) :: noco
     CHARACTER(len=1000) :: buf

     buf=ADJUSTL(line(5:len_TRIM(line)-1)
          
     READ(line,*,iostat=ios) noco%theta,noco%phi
     noco%l_soc=.TRUE.
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
   END SUBROUTINE process_soc

   SUBROUTINE process_shift(line,atompos)
     CHARACTER(len=*),INTENT(in)::line
     REAL,INTENT(INOUT) :: atompos(:,:)
     CHARACTER(len=1000) :: buf
     REAL :: shift(3)
     INTEGER :: ios,n
     
     buf=ADJUSTL(line(7:len_TRIM(line)-1)    
     READ(line,*,iostat=ios) shift
     
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
     DO n=1,SIZE(atompos,2)
        atompos(:,n)=atompos(:,)+shift
     ENDDO
   END SUBROUTINE process_shift
   SUBROUTINE process_factor(line,atompos)
     CHARACTER(len=*),INTENT(in)::line
     REAL,INTENT(INOUT) :: atompos(:,:)
     CHARACTER(len=1000) :: buf
     REAL :: factor(3)
     INTEGER :: ios,n
     
     buf=ADJUSTL(line(8:len_TRIM(line)-1)    
     READ(line,*,iostat=ios) factor
     
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
     DO n=1,SIZE(atompos,2)
        atompos(:,n)=atompos(:,)/factor
     ENDDO
   END SUBROUTINE process_factor

   SUBROUTINE process_exco(line,xcpot)
     CHARACTER(len=*),INTENT(in)::line
     TYPE(t_xcpot),INTENT(INOUT) :: xcpot
     LOGICAL::relxc
     CHARACTER(len=4) :: xctyp
     NAMELIST /exco/   xctyp, relxc 
     INTEGER :: ios
     
     READ(line,exco,iostat=ios) 
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))

     call xcpot???
   END SUBROUTINE process_exco

   SUBROUTINE process_comp(line,jspins,frcor,ctail,kcrel,gmax,gmaxxc,rkmax)
     CHARACTER(len=*),INTENT(in)::line
     INTEGER,INTENT(inout):: jspins,frcor,ctail,kcrel
     REAL,intent(inout)   :: gmax,gmaxxc,rkmax

     INTEGER :: ios
     NAMELIST /comp/   jspins, frcor, ctail, kcrel, gmax, gmaxxc, kmax
     
     READ(line,comp,iostat=ios) 
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
   END SUBROUTINE process_comp

   

      SUBROUTINE normalize_file(infh,outfh)
    !***********************************************************************
    !     reads in the file from infh
    ! and:
    !  - deletes comments
    !  - deletes empty line
    !  - combines multiple line namelists into single line
    !
    !  then the input is written to outfh
    ! 
    !***********************************************************************
    
      IMPLICIT NONE

      INTEGER, INTENT (IN)    :: infh            ! input filehandle (5)
      INTEGER, INTENT (IN)    :: outh            ! Output filehandle

      INTEGER               :: n,ios
      LOGICAL               :: building, complete
      CHARACTER(len=1000)   :: line,buffer

      
      
      !---> initialize some variables
      building = .false.
      complete = .false.

      loop: DO
         READ (infh,'(a)',IOSTAT=ios) line
         IF (ios.NE.0) EXIT !done
         
         LINE = ADJUSTL(line)
         n = SCAN(line,'!')                 ! remove end of line comments
         IF ( n>0 ) THEN
            line = line(1:n-1)
         ENDIF
         n = LEN_TRIM( line )               ! length of line without trailing blanks
         IF ( n == 0 ) CYCLE loop

         IF ( line(1:1)=='&' ) THEN         ! check if beginning of namelist
            IF (building) CALL juDFT_error ("missing end of namelist marker / in  or before line")
            building = .TRUE.
            buffer = line
            IF( line(n:n)=='/' ) complete = .TRUE.
         ELSEIF ( line(n:n)=='/' ) THEN     ! check if end of namelist
            IF (building) THEN
               complete = .TRUE.
               buffer = trim(buffer)//' '//line
            ELSE
               CALL juDFT_error ("out of place end of namelist marker / in line")
            ENDIF
         ELSEIF ( building ) THEN           ! add line to buffer
            buffer = trim(buffer)//' '//line
         ELSEIF ( n > 0 ) THEN              ! check for non empty lines outside of namelists
            buffer = line
            complete = .TRUE.
         ENDIF

         IF ( complete ) THEN
            WRITE(outfh,"(a)") TRIM(buffer)
            buffer=''
            building=.FALSE.
            complete=.FALSE.
         END IF
      END DO loop

    END SUBROUTINE normalize_file

  
    
  END MODULE m_read_inpgen_input
