MODULE m_atompar
  USE m_judft
  IMPLICIT NONE
  type t_atompar
     integer :: id = -1
     integer :: nucnumber = 0
     real    :: rmt = 0.0
     real    :: rmt_min=99.0
     integer :: jri = 0
     REAL    :: dx = 0.0
     REAL    :: bmu = -9999.0
     integer :: lmax = 0
     integer :: lnonsph = 0
     character(len=100)::lo=""
     character(len=100)::econfig=""
     character(len=100)::desc=""
   contains
     procedure :: add_defaults
  end type t_atompar

  type(t_atompar),allocatable :: atompar_list(:)
  integer :: no_of_atompars
  
contains
  subroutine add_defaults(ap)
    class(t_atompar),intent(inout)::ap

    TYPE(t_atompar):: ap_d
    INTEGER        :: n

    
    if (ap%rmt>0) then
       ap_d=find_atompar(ap%nucnumber,ap%rmt)
    else
       call judft_error("Defaults ...")
    endif

    if (ap%jri==0) ap%jri=ap_d%jri
    if (ap%dx==0) ap%jri=ap_d%dx
    if (ap%lmax==0) ap%jri=ap_d%lmax
    if (ap%lnonsph==0) ap%jri=ap_d%lnonsph
    if (ap%lo=="") ap%lo=ap_d%lo
    IF (ap%econfig=="") ap%econfig=ap_d%econfig
    
    !now generate defaults for missing values
    if(ap%jri==0) ap%jri = NINT(NINT(330*ap%rmt)*0.5)*2+1
    if(ap%dx==0) ap%dx=LOG(3200*ap%nucnumber*ap%rmt)/(ap%jri-1)
    if(ap%lmax==0) then
       if (ap%rmt<1.8) then
          ap%lmax=6
       else if(ap%rmt<2.4) then
          ap%lmax=8
       else
          ap%lmax=10
       endif
    endif
    IF(ap%lnonsph==0) ap%lnonsph=MIN( MAX( (ap%lmax-2),3 ), 8 )
    IF (ap%lnonsph>ap%lmax) THEN
       WRITE(*,*) "lnonsph had to be reduced for some atom"
       ap%lnonsph=ap%lmax
    ENDIF
    !check of magnetism
    IF (ap%bmu<-99.0) THEN
       IF (ap_d%bmu>-99.0) THEN
          ap%bmu=ap_d%bmu
       ELSE
          ap%bmu = 0.0
          IF (ap%nucnumber.EQ.24) ap%bmu = 1.0  ! Cr - Ni
          IF (ap%nucnumber.EQ.25) ap%bmu = 3.5
          IF (ap%nucnumber.EQ.26) ap%bmu = 2.2
          IF (ap%nucnumber.EQ.27) ap%bmu = 1.6
          IF (ap%nucnumber.EQ.28) ap%bmu = 1.1
          IF (ap%nucnumber.EQ.59) ap%bmu = 2.1  ! Pr - Tm
          IF (ap%nucnumber.EQ.60) ap%bmu = 3.1
          IF (ap%nucnumber.EQ.61) ap%bmu = 4.1
          IF (ap%nucnumber.EQ.62) ap%bmu = 5.1
          IF (ap%nucnumber.EQ.63) ap%bmu = 7.1
          IF (ap%nucnumber.EQ.64) ap%bmu = 7.1 
          IF (ap%nucnumber.EQ.65) ap%bmu = 6.1
          IF (ap%nucnumber.EQ.66) ap%bmu = 5.1
          IF (ap%nucnumber.EQ.67) ap%bmu = 4.1
          IF (ap%nucnumber.EQ.68) ap%bmu = 3.1
          IF (ap%nucnumber.EQ.69) ap%bmu = 2.1
       END IF
    END IF

    !make sure there are no blanks in lo
    DO n=1,len_TRIM(ap%lo)
       IF (ap%lo(n:n)==' ') THEN
          ap%lo(n:LEN(ap%lo)-n)=ap%lo(n+1:)
          ap%lo(LEN(ap%lo):LEN(ap%lo))=' '
       END IF
    ENDDO
    
    
    
  end subroutine add_defaults
    
       
  subroutine add_atompar(ap)
    TYPE(t_atompar),INTENT(in),OPTIONAL::ap
    type(t_atompar),allocatable:: tmp_list(:)

    
    
    IF (.NOT.ALLOCATED(atompar_list)) THEN
       no_of_atompars=0
       ALLOCATE(atompar_list(100))
       !Try to read default parameter files
       CALL read_params("default.econfig")
       CALL read_params("fleur.econfig")
       call read_params("my.econfig")
    else
       if (no_of_atompars==size(atompar_list)) then
          !extend the list
          call move_alloc(atompar_list,tmp_list)
          allocate(atompar_list(no_of_atompars+50))
          atompar_list(:no_of_atompars)=tmp_list
          deallocate(tmp_list)
       endif
    END IF
    IF (PRESENT(ap)) THEN
       no_of_atompars=no_of_atompars+1
       atompar_list(no_of_atompars)=ap
    ENDIF
  end subroutine add_atompar

  function find_atompar(nucnumber,rmt_max,id)result(ap)
    integer,intent(in)          :: nucnumber
    real,intent(in)             :: rmt_max
    integer,intent(in),optional :: id
    type(t_atompar)    :: ap

    integer :: n

    call add_atompar() !Make sure we have at least the defaults
    !check if there is an id given
    if (present(id)) then
       DO n=no_of_atompars,1,-1
          ap=atompar_list(n)
          IF (ap%nucnumber==nucnumber.AND.ap%id==id) THEN
             RETURN
          ENDIF
       end DO
       call judft_error("You specified a specific id for an atom but never defined that")
    end if

    !Else we search if the MT has been given for this atom
    DO n=no_of_atompars,1,-1
       ap=atompar_list(n)
       if (ap%nucnumber==nucnumber) then
          IF (ap%rmt>0.AND.ap%rmt<rmt_max) THEN
             RETURN
          ENDIF 
       endif
    enddo

    !Else we check if there is an atom definition available
    DO n=no_of_atompars,1,-1
       ap=atompar_list(n)
       if (ap%nucnumber==nucnumber) then
          IF (ap%rmt_min>5.0.OR.ap%rmt_min<=rmt_max) THEN
             ap%rmt=rmt_max
             return
          endif
       endif
    enddo
    
    call judft_error("No possible atomic parameter-set found")
  end function find_atompar


  SUBROUTINE read_params(filename)
    CHARACTER(len=*),INTENT(in)::filename

    CHARACTER(len=5)::str
    INTEGER :: id,z,jri,lmax,lnonsph,io_stat
    REAL    :: rmt,rmt_min,dx,bmu
    CHARACTER(len=100)::desc,lo,econfig
    LOGICAL :: l_exist
    TYPE(t_atompar)::ap
    NAMELIST /atom/desc,id,z,jri,lmax,lnonsph,rmt,rmt_min,dx,lo,econfig,bmu

    INQUIRE(file=filename,exist=l_exist)
    IF (.NOT.l_exist) RETURN
    
    OPEN(99,file=filename)
    DO
       READ(99,*,err=100,END=100) str
       IF (str=="&atom") THEN
          BACKSPACE(99)
          id=-1;z=0;rmt=0.0;rmt_min=99.0;jri=0;dx=0.0;lmax=0;bmu=-9999.0;lnonsph=0;lo='';econfig=''
          READ(99,atom,iostat=io_stat)
          IF (io_stat==0) THEN
             ap=t_atompar(id=id,nucnumber=z,rmt=rmt,rmt_min=rmt_min,jri=jri,dx=dx,bmu=bmu,lmax=lmax,lnonsph=lnonsph,lo=lo,econfig=econfig,desc=desc)
          ELSE
             !try old namelist
             CALL read_atom_params_old(99,ap)
          END IF
          CALL add_atompar(ap)
       ENDIF
    ENDDO
100 CLOSE(99)
  END SUBROUTINE read_params

  SUBROUTINE read_atom_params_old(fh,ap)
    !Try to read old namelist
    integer,intent(in)::fh
    TYPE(t_atompar),INTENT(out)::ap

    REAL:: id,z,rmt,dx,bmu
    INTEGER:: jri,lmax,lnonsph,ncst,nc,io_stat,nz
    CHARACTER(len=100)::econfig,lo,element,name
    
    NAMELIST /atom/ id,z,rmt,dx,jri,lmax,lnonsph,ncst,econfig,bmu,lo,element,name

    id=-9999.9;z=-1.0;rmt=0.0;dx=0.0;jri=0;lmax=0;lnonsph=0;ncst=-1;lo='';econfig='';name='';bmu=-9999.0
    
    BACKSPACE(fh)
    READ(fh,atom,iostat=io_stat)
    IF(io_stat.NE.0) THEN
       BACKSPACE(fh)
       READ(fh,*) name
       WRITE(*,*) name
       CALL judft_error("Found a &atom namelist in input that was incorrect")
    END IF
    
    !determine nz and id
    nz=-1
    IF (element.NE."") THEN
       nz=element_to_z(element)
       IF (z>-1.AND.nz.NE.FLOOR(z)) CALL judft_error("z and z of specified element differ")
    ELSE
       nz=FLOOR(z)
    ENDIF

    IF (id.NE.-9999.9) THEN
       IF (nz==-1) nz=FLOOR(id)
       id=(id-nz)*100
       IF (id>=100.OR.id<0) CALL judft_error("ID and element (or nuclear number do not fit")
    ELSE
       id=-1.0
    ENDIF
    IF (nz==-1) CALL judft_error("Neither z nor element specified")

    IF (ncst>-1) CALL judft_warn("ncst is no longer supported as input")

    IF (LEN_TRIM(econfig)==0)THEN
       ap=find_atompar(nz,rmt)
       econfig=ap%econfig
       IF (LEN_TRIM(lo)==0) lo=ap%lo
    END IF
        
    ap=t_atompar(id=INT(id),nucnumber=nz,rmt=rmt,dx=dx,jri=jri,lmax=lmax,lnonsph=lnonsph,lo=lo,bmu=bmu,econfig=econfig,desc=name)
    
    
  CONTAINS
    INTEGER FUNCTION element_to_z(element)
      USE m_constants,ONLY: namat_const
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(in)::  element
      
      CHARACTER(len=2)   :: ele,nam
      INTEGER,parameter  :: adiff=abs(ICHAR('A')-ICHAR('a'))
      INTEGER            :: n,diff
      
      ele=ADJUSTL(element)
      
      element_to_z = -1
      DO n = 0, SIZE(namat_const)-1
         nam=ADJUSTL(namat_const(n))
         diff=ABS(ICHAR(ele(1:1))-ICHAR(nam(1:1)))
         IF (diff==0.OR.diff==adiff)THEN
            diff=ABS(ICHAR(ele(2:2))-ICHAR(nam(2:2)))
            IF (diff==0.OR.diff==adiff) THEN
               element_to_z = n
               EXIT
            ENDIF
         END IF
      END DO
      
    END FUNCTION element_to_z
  END SUBROUTINE read_atom_params_old

  SUBROUTINE dump_list()
    INTEGER::n

    INTEGER :: id,z,jri,lmax,lnonsph
    REAL    :: rmt,rmt_min,dx
    CHARACTER(len=100)::desc,lo,econfig
    type(t_atompar)::ap
    NAMELIST /atom/desc,id,z,jri,lmax,lnonsph,rmt,rmt_min,dx,lo,econfig

    WRITE(6,*) "List of defined atomic parameters"
    DO n=1,no_of_atompars
       ap=atompar_list(n)
       id=ap%id
       z=ap%nucnumber
       jri=ap%jri
       lmax=ap%lmax
       lnonsph=ap%lnonsph
       rmt=ap%rmt
       rmt_min=ap%rmt_min
       dx=ap%dx
       lo=ap%lo
       econfig=ap%econfig
       desc=ap%desc
       WRITE(6,atom)
    ENDDO
    
  END SUBROUTINE dump_list
END MODULE m_atompar

