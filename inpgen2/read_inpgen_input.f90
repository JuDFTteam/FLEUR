!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_read_inpgen_input
  USE m_judft
  USE m_calculator
  IMPLICIT NONE
  PRIVATE
  PUBLIC read_inpgen_input, peekInpgenInput
CONTAINS

  SUBROUTINE read_inpgen_input(profile,atom_pos,atom_id,atom_label,kpts_str,kptsName,kptsPath,kptsBZintegration,&
                               kptsGamma,input,sym,noco,vacuum,stars,xcpot,cell,hybinp)
    !Subroutine reads the old-style input for inpgen
    USE m_atompar
    USE m_types_input
    USE m_types_sym
    USE m_types_noco
    USE m_types_vacuum
    USE m_types_stars
    USE m_types_xcpot_inbuild_nofunction
    USE m_types_cell
    USE m_types_hybinp
    USE m_constants
    USE m_process_lattice_namelist
    USE m_inv3
    USE m_types_profile

    TYPE(t_profile),INTENT(IN)     :: profile
    REAL,    ALLOCATABLE,INTENT(OUT) :: atom_pos(:, :),atom_id(:)
    CHARACTER(len=20), ALLOCATABLE,INTENT(OUT) :: atom_Label(:)
    CHARACTER(len=40),INTENT(OUT)  :: kpts_str(:)
    CHARACTER(len=40),INTENT(out)  :: kptsName(:)
    CHARACTER(len=500),INTENT(out) :: kptsPath(:)
    INTEGER,INTENT(OUT)            :: kptsBZintegration(:)
    LOGICAL,INTENT(OUT)            :: kptsGamma(:)
    TYPE(t_input),INTENT(out)      :: input
    TYPE(t_sym),INTENT(OUT)        :: sym
    TYPE(t_noco),INTENT(OUT)       :: noco
    TYPE(t_vacuum),INTENT(OUT)     :: vacuum
    TYPE(t_stars),INTENT(OUT)      :: stars
    TYPE(t_xcpot_inbuild_nf),INTENT(OUT)    :: xcpot
    TYPE(t_cell),INTENT(OUT)       :: cell
    TYPE(t_hybinp),INTENT(OUT)     :: hybinp



    !locals
    REAL                :: a1(3),a2(3),a3(3),aa,SCALE(3),mat(3,3),cart_mat(3,3),det,temp
    INTEGER             :: ios,n,i, iKpts
    CHARACTER(len=100)  :: filename
    LOGICAL             :: l_exist
    LOGICAL             :: cartesian=.false.
    CHARACTER(len=16384):: line
    TYPE(t_atompar)     :: ap

    iKpts = 0

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
    CLOSE(97)

    REWIND(98)
    READ(98,"(a)",iostat=ios) input%comment
    IF (input%comment(1:1)=='&') THEN
       input%comment="No title"
       BACKSPACE(98)
    ENDIF

    IF(TRIM(ADJUSTL(profile%profileName)).NE."default") THEN
       input%tkb = profile%fermiSmearing
    END IF

    aa=0.0
    input%jspins=0
    readloop: DO WHILE(ios==0)
       READ(98,"(a)",iostat=ios) line
       IF (ios.NE.0) EXIT
       IF (line(1:1)=="&") THEN
          !process the namelist
          SELECT CASE(line(2:5)) !e.g. atom
          CASE ('latt')
             CALL process_lattice(line,a1,a2,a3,aa,scale,mat,cart_mat)
          CASE('inpu')
             CALL process_input(line,input%film,sym%symor,cartesian,hybinp%l_hybrid)
          CASE('atom')
             CALL read_atom_params_old(98,ap,profile)
             CALL add_atompar(ap)
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
          CASE('expe')
             CALL process_expert(line,input%gw,cell%primCellZ)
          CASE('kpt ')
             iKpts = iKpts + 1
             CALL process_kpts(line,kpts_str(iKpts),kptsName(iKpts),kptsPath(iKpts),kptsBZintegration(iKpts),kptsGamma(ikpts),input%tkb)
             IF(TRIM(ADJUSTL(kptsName(iKpts))).EQ.'') THEN
                IF(TRIM(ADJUSTL(kptsPath(iKpts))).EQ.'') THEN
                   WRITE(kptsName(iKpts),'(a,i0)') "default-", iKpts
                ELSE
                   WRITE(kptsName(iKpts),'(a,i0)') "path-", iKpts
                END IF
             END IF
          CASE('film')
             CALL process_film(line,vacuum%dvac,cell%amat(3,3))
          CASE('gen ','sym ')
             CALL judft_error("Specifying the symmetries no longer supported in inpgen")
          CASE('end')
             EXIT readloop
          CASE default
             CALL judft_error(("Unknown input in:"//trim(line)))
          END SELECT
       ELSE
          IF (aa.NE.0) THEN
             !cell was set already, so list of atoms follow
             if (allocated(atom_pos)) call judft_error("Input error: "//TRIM(line))
             READ(line,*,iostat=ios) n
             IF (ios.NE.0) CALL judft_error(("Surprising error in reading input: "//trim(line)))
             ALLOCATE(atom_pos(3,n),atom_label(n),atom_id(n))
             DO i=1,n
                READ(98,"(a)",iostat=ios) line
                IF (ios.NE.0) CALL judft_error(("List of atoms not complete: "//trim(line)))
                atom_id(i)=evaluatefirst(line)
                atom_pos(1,i)=evaluatefirst(line)
                atom_pos(2,i)=evaluatefirst(line)
                atom_pos(3,i)=evaluatefirst(line)
                IF(TRIM(ADJUSTL(line)).NE.'') THEN
                   atom_Label(i) = TRIM(ADJUSTL(line))
                ELSE
                   WRITE(atom_Label(i),'(i0)') i
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
                vacuum%dvac=0.01
                !vacuum%dvac = ABS(a3(3)) ! This is later set to the real default by the chkmt result
             END IF
             READ(98,"(a)",iostat=ios) line
             IF (ios.NE.0) CALL judft_error(("Error reading bravais matrix"))
             aa=evaluatefirst(line)
             READ(98,"(a)",iostat=ios) line
             IF (ios.NE.0) CALL judft_error(("Error reading bravais matrix"))
             SCALE(1)=evaluatefirst(line);SCALE(2)=evaluatefirst(line);SCALE(3)=evaluatefirst(line)
             mat=0.0
          ENDIF
       ENDIF
    END DO readloop
    
    IF ((ikpts.EQ.0).AND.(profile%kPDen.GT.0.0)) THEN
       iKpts = iKpts + 1
       line = ''
       WRITE(line,'(a,f12.4,a,f15.10,a)') '&kpt gamma=T den=', profile%kPDen, ' tkb=', profile%fermiSmearing, ' /'
       WRITE(*,*) 'kPDen: ', profile%kPDen
       WRITE(*,'(a)') TRIM(line)
       CALL process_kpts(line,kpts_str(iKpts),kptsName(iKpts),kptsPath(iKpts),kptsBZintegration(iKpts),kptsGamma(ikpts),input%tkb)
    END IF

    IF (.NOT.ALLOCATED(atom_pos).OR.SUM(ABS(a1))==0.0) CALL judft_error("input not complete")

    !transform hex->trig
    IF (ABS(mat(1,1)).GT.0.0000001) THEN
        !unscaled matrices... (scaled setup later)
        cell%amat(:,1) = a1(:)
        cell%amat(:,2) = a2(:)
        cell%amat(:,3) = a3(:)
        CALL inv3(cell%amat,cell%bmat,det)
        DO n = 1, SIZE(atom_pos,2)
          atom_pos(:,n) = MATMUL(cell%bmat,MATMUL(mat,atom_pos(:,n)))
        ENDDO
    ENDIF
    !Transform in case of scaled cartesian input
    IF (cartesian) THEN
      if (all(abs(cart_mat)<0.01)) call judft_error("Cartesian='t' not possible for your lattice")
      CALL inv3(cart_mat,cell%bmat,det)
      DO n = 1, SIZE(atom_pos,2)
          atom_pos(:,n) = MATMUL(cell%bmat,atom_pos(:,n))
       ENDDO
    ENDIF
    DO i = 1, 3
       IF (SCALE(i).LT.0.0) SCALE(i) = SQRT(-SCALE(i))
    END DO

    !set the cell
    cell%amat(:,1) = aa*SCALE(:)*a1(:)
    cell%amat(:,2) = aa*SCALE(:)*a2(:)
    cell%amat(:,3) = aa*SCALE(:)*a3(:)

    !convert to right-handed unit cell if it is left-handed so far
    CALL inv3(cell%amat,cell%bmat,det)
    IF(det.LT.0.0) THEN
       cell%amat(:,1) = aa*SCALE(:)*a2(:)
       cell%amat(:,2) = aa*SCALE(:)*a1(:)
       DO n = 1, SIZE(atom_pos,2)
          temp = atom_pos(1,n)
          atom_pos(1,n) = atom_pos(2,n)
          atom_pos(2,n) = temp
       END DO
       WRITE(*,*) 'Provided unit cell is left-handed. Converting it to right-handed system by exchanging a1 and a2'
       WRITE(oUnit,*) 'Provided unit cell is left-handed. Converting it to right-handed system by exchanging a1 and a2'
    END IF

    CALL cell%init(0.0)

    CLOSE(98)

  END SUBROUTINE read_inpgen_input

  SUBROUTINE peekInpgenInput(numKpts,numKptsPath)

    INTEGER, INTENT(INOUT) :: numKpts
    INTEGER, INTENT(INOUT) :: numKptsPath

    INTEGER             :: ios
    LOGICAL             :: l_exist,gamma
    CHARACTER(len=100)  :: filename
    CHARACTER(len=16384):: line

    CHARACTER(len=40)  :: kpts_str
    CHARACTER(len=40)  :: kptsName
    CHARACTER(len=500) :: kptsPath
    INTEGER            :: bz_integration_out
    REAL               :: tkb

    numKpts = 0
    numKptsPath = 0
    ios = 0

    filename=juDFT_string_for_argument("-f")
    INQUIRE(file=filename,exist=l_exist)
    IF (.NOT.l_exist) CALL judft_error("Input file specified is not readable")
    OPEN(97,file=filename)
    !OPEN(98,status='scratch')
    OPEN(98,file="scratch")
    CALL normalize_file(97,98)
    CLOSE(97)

    REWIND(98)

    DO WHILE(ios==0)
       READ(98,"(a)",iostat=ios) line
       IF (ios.NE.0) EXIT
       IF (line(1:1)=="&") THEN
          SELECT CASE(line(2:5)) !e.g. atom
          CASE('kpt ')
             numKpts = numKpts + 1
             CALL process_kpts(line,kpts_str,kptsName,kptsPath,bz_integration_out,gamma,tkb)
             IF(kptsPath.NE.'') numKptsPath = numKptsPath + 1
          END SELECT
       END IF
    END DO

    CLOSE(98)

  END SUBROUTINE peekInpgenInput


  SUBROUTINE process_kpts(line,kpts_str,kptsName,kptsPath,bz_integration_out,kptsGamma,tkb)
    USE m_constants
    CHARACTER(len=*),INTENT(in)::line
    CHARACTER(len=40),INTENT(out)::kpts_str
    CHARACTER(len=40),INTENT(out)::kptsName
    CHARACTER(len=500),INTENT(out)::kptsPath
    INTEGER,INTENT(inout)::bz_integration_out
    LOGICAL,INTENT(out)::kptsGamma
    REAL,INTENT(inout):: tkb

    LOGICAL :: tria, gamma
    INTEGER :: div1,div2,div3,nkpt, numSpecifications
    CHARACTER(len=5) :: bz_integration
    CHARACTER(len=40) :: name
    CHARACTER(len=500) :: path
    REAL    :: den
    NAMELIST /kpt/nkpt,div1,div2,div3,tkb,bz_integration,gamma,tria,den,path,name
    div1=0;div2=0;div3=0;nkpt=0;den=0.0
    bz_integration = 'hist'
    name = ''
    path = ''
    tria=.FALSE.
    gamma=.FALSE.
    READ(line,kpt)
    kpts_str=''

    numSpecifications = 0
    IF (den.GT.0.0) numSpecifications = numSpecifications + 1
    IF (nkpt.NE.0) numSpecifications = numSpecifications + 1
    IF (ALL([div1,div2,div3]>0)) numSpecifications = numSpecifications + 1
    IF (numSpecifications.GT.1) THEN
       WRITE(*,'(a)') "Ambiguous specification of k-point set:"
       WRITE(*,'(a)') TRIM(line)
       CALL juDFT_error("Ambiguous specification of k-point set.", calledby="read_inpgen_input",&
                        hint="Use only one of nkpt, den, (div1,div2,div3)!")
    END IF

    IF (den>0.0) THEN
       WRITE(kpts_str,"(a,f0.6)") "den=",den
    ELSEIF((nkpt>0).AND.(path.EQ.'')) THEN
       WRITE(kpts_str,"(a,i0)") "nk=",nkpt
    ELSEIF((nkpt>0).AND.(path.NE.'')) THEN
       WRITE(kpts_str,"(a,i0)") "band=",nkpt
    ELSEIF(ALL([div1,div2,div3]>0)) THEN
       WRITE(kpts_str,"(a,i0,a,i0,a,i0)") "grid=",div1,",",div2,",",div3
    END IF
    SELECT CASE(TRIM(ADJUSTL(bz_integration)))
    CASE('hist')
       bz_integration_out = BZINT_METHOD_HIST
    CASE('gauss')
       bz_integration_out = BZINT_METHOD_GAUSS
    CASE('tria')
       bz_integration_out = BZINT_METHOD_TRIA
    CASE('tetra')
       bz_integration_out = BZINT_METHOD_TETRA
    CASE DEFAULT
       CALL judft_error("No valid bz_integration mode",calledby="process_kpts")
    END SELECT
    IF(tria.AND.bz_integration_out==BZINT_METHOD_HIST) THEN
       bz_integration_out = BZINT_METHOD_TRIA
    ELSE IF(tria.AND.bz_integration_out /= BZINT_METHOD_HIST) THEN
       CALL juDFT_warn("You specified both tria and bz_integration in the input,\\&
                       tria will be ignored",calledby="process_kpts")
    ENDIF
    kptsName = name
    kptsPath = path
    kptsGamma = gamma

  END SUBROUTINE process_kpts

  SUBROUTINE process_input(line,film,symor,cartesian,hybinp)
    CHARACTER(len=*),INTENT(in)::line
    LOGICAL,INTENT(out)::film,symor,hybinp


    INTEGER :: ios
    LOGICAL :: cartesian, cal_symm, checkinp,inistop,oldfleur
    NAMELIST /input/ film, cartesian, cal_symm, checkinp, inistop,&
         symor, oldfleur, hybinp
    cartesian=.FALSE.
    cal_symm=.FALSE.
    oldfleur=.FALSE.
    checkinp=.FALSE.
    READ(line,input,iostat=ios)
    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))
    IF (ANY([cal_symm, checkinp,oldfleur])) CALL judft_error("Switches cal_symm, checkinp,oldfleur no longer supported")
  END SUBROUTINE process_input

  SUBROUTINE process_qss(line,noco)
    USE m_types_noco
    CHARACTER(len=*),INTENT(in)::line
    TYPE(t_noco),INTENT(INOUT) :: noco
    CHARACTER(len=1000) :: buf
    INTEGER :: ios

    buf=ADJUSTL(line(5:LEN_TRIM(line)-1))

    READ(buf,*,iostat=ios) noco%qss_inp
    noco%l_ss=.TRUE.
    noco%l_noco=.TRUE.
    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))
  END SUBROUTINE process_qss

  SUBROUTINE process_soc(line,noco)
    USE m_types_noco
    CHARACTER(len=*),INTENT(in)::line
    TYPE(t_noco),INTENT(INOUT) :: noco
    CHARACTER(len=1000) :: buf
    INTEGER :: ios

    buf=ADJUSTL(line(5:LEN_TRIM(line)-1))

    READ(buf,*,iostat=ios) noco%theta_inp,noco%phi_inp
    noco%l_soc=.TRUE.
    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))
  END SUBROUTINE process_soc

  SUBROUTINE process_film(line,dvac,dtild)
    CHARACTER(len=*),INTENT(in):: line
    REAL,INTENT(out)           :: dvac,dtild
    INTEGER :: ios
    NAMELIST /film/   dvac, dtild
    READ(line,film,iostat=ios)
    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))
  END SUBROUTINE process_film

  SUBROUTINE process_shift(line,atompos)
    CHARACTER(len=*),INTENT(in)::line
    REAL,INTENT(INOUT) :: atompos(:,:)
    CHARACTER(len=1000) :: buf
    REAL :: shift(3)
    INTEGER :: ios,n

    buf=ADJUSTL(line(7:LEN_TRIM(line)-1))
    READ(buf,*,iostat=ios) shift

    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))
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

    buf=ADJUSTL(line(8:LEN_TRIM(line)-1))
    READ(buf,*,iostat=ios) factor

    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))
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

    relxc=.FALSE.
    xctyp='pbe'
    READ(line,exco,iostat=ios)
    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))

    xcpot%l_inbuild=.TRUE.
    xcpot%inbuild_name=xctyp
    xcpot%l_relativistic=relxc
    CALL xcpot%init(1) !Is it OK to use ntype=1 here??
  END SUBROUTINE process_exco

  SUBROUTINE process_comp(line,jspins,frcor,ctail,kcrel,gmax,gmaxxc,kmax)
    CHARACTER(len=*),INTENT(in)::line
    INTEGER,INTENT(inout):: jspins,kcrel
    LOGICAL,INTENT(inout):: frcor,ctail
    REAL,INTENT(inout)   :: gmax,gmaxxc,kmax

    INTEGER :: ios
    NAMELIST /comp/   jspins, frcor, ctail, kcrel, gmax, gmaxxc, kmax

    READ(line,comp,iostat=ios)
    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))
  END SUBROUTINE process_comp

  SUBROUTINE process_expert(line,gw,primCellZ)
    USE m_types_xcpot_inbuild_nofunction
    CHARACTER(len=*),INTENT(in)::line
    INTEGER, INTENT(INOUT) :: gw
    REAL, INTENT(OUT) :: primCellZ
    INTEGER :: spex
    INTEGER :: ios
    NAMELIST /expert/ spex, primCellZ
    primCellZ = 0.0


    spex = 0
    READ(line,expert,iostat=ios)
    IF (ios.NE.0) CALL judft_error(("Error reading:" //TRIM(line)))

    gw = spex
  END SUBROUTINE process_expert

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
    building = .FALSE.
    complete = .FALSE.

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
             buffer = TRIM(buffer)//' '//line
          ELSE
             CALL juDFT_error ("out of place end of namelist marker / in line")
          ENDIF
       ELSEIF ( building ) THEN           ! add line to buffer
          buffer = TRIM(buffer)//' '//line
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
    IF(LEN_TRIM(buffer).NE.0) THEN
       WRITE(outfh,'(a)') TRIM(buffer)
    END IF

  END SUBROUTINE normalize_file



END MODULE m_read_inpgen_input
