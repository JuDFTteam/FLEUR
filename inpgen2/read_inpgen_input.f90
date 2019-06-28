!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_read_inpgen_input
  USE m_judft
  USE m_calculator
  IMPLICIT NONE
  private
  public read_inpgen_input
CONTAINS

  SUBROUTINE read_inpgen_input(atom_pos,atom_id,atom_label,kpts_str,&
       input,sym,noco,vacuum,stars,xcpot,cell,hybrid)
    !Subroutine reads the old-style input for inpgen
    use m_atompar
    USE m_types_input
    USE m_types_sym
    USE m_types_noco
    USE m_types_vacuum
    USE m_types_stars
    USE m_types_xcpot_inbuild_nofunction
    USE m_types_cell
    USE m_types_hybrid
    USE m_process_lattice_namelist
    
    use m_inv3
    REAL,    ALLOCATABLE,INTENT(OUT) :: atom_pos(:, :),atom_id(:)
    CHARACTER(len=20), ALLOCATABLE,INTENT(OUT) :: atom_Label(:)
    character(len=40),INTENT(OUT):: kpts_str
    TYPE(t_input),INTENT(out)    :: input
    TYPE(t_sym),INTENT(OUT)      :: sym
    TYPE(t_noco),INTENT(OUT)     :: noco
    TYPE(t_vacuum),INTENT(OUT)   :: vacuum
    TYPE(t_stars),INTENT(OUT)    :: stars
    TYPE(t_xcpot_inbuild_nf),INTENT(OUT)    :: xcpot
    TYPE(t_cell),INTENT(OUT)     :: cell
    TYPE(t_hybrid),INTENT(OUT)   :: hybrid
    
    
    !locals
    REAL                :: a1(3),a2(3),a3(3),aa,SCALE(3),mat(3,3),det
    integer             :: ios,n,i
    CHARACTER(len=100)  :: filename
    logical             :: l_exist
    CHARACTER(len=16384):: line
    type(t_atompar)     :: ap
     
    !
    !read file and create normalized scratch file
    !
    filename=juDFT_string_for_argument("-f")
    INQUIRE(file=filename,exist=l_exist)
    IF (.NOT.l_exist) CALL judft_error("Input file specified is not readable")
    OPEN(97,file=filename)
    !OPEN(98,status='scratch')
    OPEN(98,file="scratch")
    CALL  normalize_file(97,98)
    close(97)
    
    REWIND(98)
    READ(98,"(a)",iostat=ios) input%comment
    if (input%comment(1:1)=='&') then
       input%comment="No title"
       backspace(98)
    endif

    aa=0.0
    DO WHILE(ios==0)
       READ(98,"(a)",iostat=ios) line
       IF (ios.NE.0) EXIT
       IF (line(1:1)=="&") THEN
          !process the namelist
          SELECT CASE(line(2:5)) !e.g. atom
          CASE ('latt')
             CALL process_lattice(line,a1,a2,a3,aa,scale,mat)
          CASE('inpu')
             CALL process_input(line,input%film,sym%symor,hybrid%l_hybrid)
          CASE('atom')
             CALL read_atom_params_old(98,ap)
             call add_atompar(ap)
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
             CALL process_kpts(line,kpts_str,input%tria,input%tkb)
          CASE('film')
             CALL process_film(line,vacuum%dvac,cell%amat(3,3))
          CASE('gen ','sym ')
             CALL judft_error("Specifying the symmetries no longer supported in inpgen")
          CASE default
             CALL judft_error(("Unkown input in:"//line))
          END SELECT
       ELSE
          IF (aa.ne.0) THEN
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
             a1(1)=evaluatefirst(line);a1(2)=evaluatefirst(line);a1(3)=evaluatefirst(line)
             READ(98,"(a)",iostat=ios) line
             IF (ios.NE.0) CALL judft_error(("Error reading bravais matrix"))
             a2(1)=evaluatefirst(line);a2(2)=evaluatefirst(line);a2(3)=evaluatefirst(line)
             READ(98,"(a)",iostat=ios) line
             IF (ios.NE.0) CALL judft_error(("Error reading bravais matrix"))
             a3(1)=evaluatefirst(line);a3(2)=evaluatefirst(line);a3(3)=evaluatefirst(line)
             vacuum%dvac=evaluatefirst(line)
             IF (input%film.AND.(vacuum%dvac <= 0.00)) THEN
                WRITE(*,*)'Film calculation but no reasonable dVac provided'
                WRITE(*,*)'Setting default for dVac'
                vacuum%dvac = ABS(a3(3)) ! This is later set to the real default by the chkmt result
             END IF
             READ(98,"(a)",iostat=ios) line
             IF (ios.NE.0) CALL judft_error(("Error reading bravais matrix"))
             aa=evaluatefirst(line)
             READ(98,"(a)",iostat=ios) line
             IF (ios.NE.0) CALL judft_error(("Error reading bravais matrix"))
             scale(1)=evaluatefirst(line);scale(2)=evaluatefirst(line);scale(3)=evaluatefirst(line)
             mat=0.0
          ENDIF
       ENDIF
    END DO

    IF (.NOT.ALLOCATED(atom_pos).OR.SUM(ABS(a1))==0.0) CALL judft_error("input not complete")
    !transform hex->trig
    IF (abs(mat(1,1)).GT.0.0000001) THEN
       cell%amat(:,1) = a1(:)
       cell%amat(:,2) = a2(:)
       cell%amat(:,3) = a3(:)
       call inv3(cell%amat,cell%bmat,det)
       DO n = 1, size(atom_pos,2)
          atom_pos(:,n) = matmul(cell%bmat,matmul(mat,atom_pos(:,n)))
       ENDDO
    ENDIF
    
    !set the cell
    cell%amat(:,1) = aa*scale(:)*a1(:)
    cell%amat(:,2) = aa*scale(:)*a2(:)
    cell%amat(:,3) = aa*scale(:)*a3(:)
    CALL cell%init()
    
  END SUBROUTINE read_inpgen_input


  SUBROUTINE process_kpts(line,kpts_str,tria,tkb)
      CHARACTER(len=*),INTENT(in)::line
      CHARACTER(len=40),INTENT(out)::kpts_str
      logical,intent(inout)::tria
      real,intent(inout):: tkb


      integer :: div1,div2,div3,nkpt
      real    :: den
      NAMELIST /kpt/nkpt,div1,div2,div3,tkb,tria,den
      div1=0;div2=0;div3=0;nkpt=0;den=0.0

      kpts_str=''
      if (den>0.0) THEN
         write(kpts_str,"(a,f0.6)") "den=",den
      elseif(nkpt>0) then
         write(kpts_str,"(a,i0)") "nk=",nkpt
      elseif(all([div1,div2,div3]>0)) then
         write(kpts_str,"(a,i0,a,i0,a,i0)") "grid=",div1,",",div2,",",div3
      end if
    end SUBROUTINE process_kpts
      
    SUBROUTINE process_input(line,film,symor,hybrid)
      CHARACTER(len=*),INTENT(in)::line
      LOGICAL,INTENT(out)::film,symor,hybrid

      INTEGER :: ios
      LOGICAL :: cartesian, cal_symm, checkinp,inistop,oldfleur
      NAMELIST /input/ film, cartesian, cal_symm, checkinp, inistop,&
                      symor, oldfleur, hybrid
      cartesian=.FALSE.
      cal_symm=.FALSE.
      oldfleur=.FALSE.
     READ(line,input,iostat=ios)
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
     IF (ANY([cal_symm, checkinp,oldfleur])) CALL judft_error("Switches cal_symm, checkinp,oldfleur no longer supported")
   END SUBROUTINE process_input
   
   SUBROUTINE process_qss(line,noco)
     USE m_types_noco
     CHARACTER(len=*),INTENT(in)::line
     TYPE(t_noco),INTENT(INOUT) :: noco
     CHARACTER(len=1000) :: buf
     INTEGER :: ios
 
     buf=ADJUSTL(line(5:len_TRIM(line)-1))
          
     READ(buf,*,iostat=ios) noco%qss
     noco%l_ss=.TRUE.
     noco%l_noco=.TRUE.
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
   END SUBROUTINE process_qss

   SUBROUTINE process_soc(line,noco)
     USE m_types_noco
     CHARACTER(len=*),INTENT(in)::line
     TYPE(t_noco),INTENT(INOUT) :: noco
     CHARACTER(len=1000) :: buf
     INTEGER :: ios
 
     buf=ADJUSTL(line(5:len_TRIM(line)-1))
          
     READ(buf,*,iostat=ios) noco%theta,noco%phi
     noco%l_soc=.TRUE.
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
   END SUBROUTINE process_soc

   SUBROUTINE process_film(line,dvac,dtild)
     CHARACTER(len=*),INTENT(in):: line
     REAL,INTENT(out)           :: dvac,dtild
     INTEGER :: ios
     NAMELIST /film/   dvac, dtild
     READ(line,film,iostat=ios)
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
   END SUBROUTINE process_film

   SUBROUTINE process_shift(line,atompos)
     CHARACTER(len=*),INTENT(in)::line
     REAL,INTENT(INOUT) :: atompos(:,:)
     CHARACTER(len=1000) :: buf
     REAL :: shift(3)
     INTEGER :: ios,n
     
     buf=ADJUSTL(line(7:len_TRIM(line)-1))    
     READ(buf,*,iostat=ios) shift
     
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
     DO n=1,SIZE(atompos,2)
        atompos(:,n)=atompos(:,n)+shift
     ENDDO
   END SUBROUTINE process_shift

   SUBROUTINE process_factor(line,atompos)
     CHARACTER(len=*),INTENT(in)::line
     REAL,INTENT(INOUT) :: atompos(:,:)
     CHARACTER(len=1000) :: buf
     REAL :: factor(3)
     INTEGER :: ios,n
     
     buf=ADJUSTL(line(8:len_TRIM(line)-1))    
     READ(buf,*,iostat=ios) factor
     
     IF (ios.NE.0) CALL judft_error(("Error reading:"//line))
     DO n=1,SIZE(atompos,2)
        atompos(:,n)=atompos(:,n)/factor
     ENDDO
   END SUBROUTINE process_factor

   SUBROUTINE process_exco(line,xcpot)
     USE m_types_xcpot_inbuild_nofunction
     CHARACTER(len=*),INTENT(in)::line
     TYPE(t_xcpot_inbuild_nf),INTENT(INOUT) :: xcpot
     LOGICAL::relxc
     CHARACTER(len=4) :: xctyp
     NAMELIST /exco/   xctyp, relxc 
     INTEGER :: ios
     
     relxc=.false.
     xctyp='pbe'
     READ(line,exco,iostat=ios) 
     IF (ios.NE.0) CALL judft_error(("Error reading:"//trim(line)))

     xcpot%l_inbuild=.true.
     xcpot%inbuild_name=xctyp
     xcpot%l_relativistic=relxc
     call xcpot%init(1) !Is it OK to use ntype=1 here??
   END SUBROUTINE process_exco

   SUBROUTINE process_comp(line,jspins,frcor,ctail,kcrel,gmax,gmaxxc,kmax)
     CHARACTER(len=*),INTENT(in)::line
     INTEGER,INTENT(inout):: jspins,kcrel
     LOGICAL,INTENT(inout):: frcor,ctail
     REAL,intent(inout)   :: gmax,gmaxxc,kmax

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
      INTEGER, INTENT (IN)    :: outfh            ! Output filehandle

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
