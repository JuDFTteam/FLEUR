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


    ALLOCATE ( mmrot(3,3,48), ttr(3,48) )
    ALLOCATE ( atompos(3,natmax),atomid(natmax) )
    ALLOCATE (atomLabel(natmax))
    atomLabel = ''

    OPEN (bfh,file='bfh.txt',form='formatted',status='unknown')
    IF (PRESENT(filename)) OPEN(5,file=filename)
    
    noco%l_ss = .FALSE.
    vacuum%dvac=0.0
    noco%l_soc=.FALSE.

    
    !default file handlers
    infh = 5
    errfh = 6 ; warnfh = 6 ; dbgfh = 6 ; outfh = 6;symfh=97

    
    nline=0
    CALL struct_input(&
                       infh,errfh,warnfh,symfh,'sym    ',bfh,&
                       natmax,48,&
                       nline,size(buffer),buffer,&
                       title,input%film,cal_symm,checkinp,sym%symor,&
                       cartesian,oldfleur,a1,a2,a3,vacuum%dvac,aa,scale,i_c,&
                       factor,natin,atomid,atompos,ngen,mmrot,ttr,atomLabel,&
                       l_hyb,noco%l_soc,noco%l_ss,noco%theta,noco%phi,noco%qss,inistop)

    !Check output
    IF (.NOT.cal_symm) CALL judft_error("Reading of symmetry no longer supported")
    IF (checkinp.OR.inistop) CALL judft_warn("checkinp and inistop no longer supported")
    IF (cartesian) CALL judft_error("Scaled Cartesian coordinates no longer supported")
    IF (l_hyb) CALL judft_warn("Hybrid option no longer supported")
    if (oldfleur.and..not.input%film)  CALL judft_warn("oldfleur only in film setups")

    !Generate amat
    amat(:,1)=aa*SCALE(:)*a1(:)
    amat(:,2)=aa*SCALE(:)*a2(:)
    amat(:,3)=aa*SCALE(:)*a3(:)
    !Generate list of atoms
    ALLOCATE(atom_pos(3,natin),atom_id(natin),atom_label(natin))
    atom_pos=atompos(:,:natin)
    atom_id=atomid(:,:natin)
    atom_label=atomlabel(:,:natin)
    !title
    DO i = 1, 10
       j = (i-1) * 8 + 1
       input%comment(i) = title(j:j+7)
    ENDDO
    IF (.NOT.input%film) vacuum%dvac=a3(3)
    dtild=0.0
    input%l_inpXML = .TRUE.
    
    CALL lapw_input(&
                     infh,nline,xl_buffer,bfh,buffer,&
                     input%jspins,input%kcrel,obsolete%ndvgrd,kpts%nkpt,div,kpts%kPointDensity,&
                     input%frcor,input%ctail,obsolete%chng,input%tria,input%rkmax,stars%gmax,xcpot%gmaxxc,&
                     vacuum%dvac,dtild,input%tkb,namex,relcor)

      CLOSE(bfh)

      !Read the &atom namelists and put into atompar as defaults
      
      CALL read_params("bfh.txt")

      OPEN(bfh,"bfh,txt")
      CLOSE(bfh,status='delete')
      
      IF (PRESENT(filename)) CLOSE(5)
      
      
    END SUBROUTINE read_inpgen_input
  END MODULE m_read_inpgen_input
