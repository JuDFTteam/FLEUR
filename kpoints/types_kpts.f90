!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts
  USE m_judft
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  TYPE,EXTENDS(t_fleurinput_base):: t_kpts
     character(len=20)     :: name="default"
     INTEGER               :: nkpt=0
     INTEGER               :: ntet=0
     LOGICAL               :: l_gamma=.FALSE.
     !(3,nkpt) k-vectors internal units
     REAL,ALLOCATABLE      :: bk(:,:)
     !(nkpts) weights
     REAL,ALLOCATABLE      :: wtkpt(:)
     INTEGER               :: nkptf=0   !<k-vectors in full BZ
     REAL   ,ALLOCATABLE   :: bkf(:,:)
     INTEGER,ALLOCATABLE   :: bkp(:)
     INTEGER,ALLOCATABLE   :: bksym(:)
     INTEGER                       :: numSpecialPoints=0
     INTEGER, ALLOCATABLE          :: specialPointIndices(:)
     CHARACTER(LEN=50),ALLOCATABLE :: specialPointNames(:)
     REAL   ,ALLOCATABLE           :: specialPoints(:,:)
     INTEGER,ALLOCATABLE           :: ntetra(:,:)
     REAL   ,ALLOCATABLE           :: voltet(:)
     REAL   ,ALLOCATABLE           :: sc_list(:,:) !list for all information about folding of bandstructure (need for unfoldBandKPTS)((k(x,y,z),K(x,y,z),m(g1,g2,g3)),(nkpt),k_original(x,y,z))
   CONTAINS
     PROCEDURE :: init
     !procedure :: read_xml
     PROCEDURE :: add_special_line
     PROCEDURE :: print_xml
     PROCEDURE :: read_xml
  ENDTYPE t_kpts

  PUBLIC :: t_kpts
CONTAINS


  SUBROUTINE read_xml(this,xml)
    USE m_types_xml
    USE m_calculator
    CLASS(t_kpts),INTENT(out):: this
    TYPE(t_xml)              :: xml
    CHARACTER(len=*)  ::name !TODO
    
    INTEGER:: number_sets,n
    CHARACTER(len=200)::str,path,path2
    


     number_sets = xml%GetNumberOfNodes('/fleurInput/calculationSetup/bzIntegration/kPointList')

     DO n=1,number_sets
        WRITE(path,"(a,i0,a)") '/fleurInput/calculationSetup/bzIntegration/kPointList[',n,']'
        IF(TRIM(ADJUSTL(name))==xml%GetAttributeValue(TRIM(path)//'/@name')) EXIT
     enddo
     IF (n>number_sets) CALL judft_error(("No kpoints named:"//TRIM(name)//" found"))
     this%nkpt=evaluateFirstOnly(xml%GetAttributeValue(TRIM(path)//'/@count'))

     this%numSpecialPoints=xml%GetNumberOfNodes(TRIM(path)//"specialPoint")

     IF (this%numSpecialPoints>0) THEN
        If (this%numSpecialPoints<2) CALL judft_error(("Giving less than two sepcial points make no sense:"//TRIM(name)))
        ALLOCATE(this%specialPoints(3,this%numSpecialPoints))
        ALLOCATE(this%specialPointNames(this%numSpecialPoints))
        DO n=1,this%numSpecialPoints
           WRITE(path2,"(a,a,i0,a)") path,"specialPoint[",n,"]"
           this%specialPointNames(n)=xml%getAttributeValue(path2//"/@name")
           str=xml%getAttributeValue(path2)
           this%specialPoints(1,n) = evaluatefirst(str)
           this%specialPoints(2,n) = evaluatefirst(str)
           this%specialPoints(3,n) = evaluatefirst(str)
        ENDDO
     ELSE
        n=xml%GetNumberOfNodes(TRIM(path)//'kPoint')
        IF (n.NE.this%nkpt) CALL judft_error(("Inconsistent number of k-points in:"//name))
        ALLOCATE(this%bk(3,this%nkpt))
        ALLOCATE(this%wtkpt(this%nkpt))
        DO n=1,this%nkpt
           WRITE(path2,"(a,a,i0,a)") path,"/kPoint[",n,"]"
           this%wtkpt(n) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(path2)//'/@weight'))
           str= xml%getAttributeValue(path2)
           this%bk(1,n) = evaluatefirst(str)
           this%bk(2,n) = evaluatefirst(str)
           this%bk(3,n) = evaluatefirst(str)
        ENDDO
        this%ntet=xml%GetNumberOfNodes(TRIM(path)//'/tetraeder/tet')
        ALLOCATE(this%voltet(this%ntet),this%ntetra(4,this%ntet))
        DO n=1,this%ntet
           WRITE(path2,"(a,a,i0,a)") path,"/tetraeder/tet[",n,"]"
           this%voltet(n)=evaluateFirstOnly(xml%GetAttributeValue(TRIM(path2)//'/@vol'))
           str= xml%getAttributeValue(path2)
           READ(str,*) this%ntetra(:,n)
        ENDDO
     END IF
   END SUBROUTINE read_xml

  SUBROUTINE print_xml(kpts,fh,filename)
    CLASS(t_kpts),INTENT(in   ):: kpts
    INTEGER,INTENT(in)         :: fh
    CHARACTER(len=*),INTENT(in),OPTIONAL::filename

    INTEGER :: n
    LOGICAL :: l_exist

    IF (PRESENT(filename)) THEN
       INQUIRE(file=filename,exist=l_exist)
       IF (l_exist) THEN
          OPEN(fh,file=filename,action="write",position="append")
       ELSE
          OPEN(fh,file=filename,action="write")
       END IF
    ENDIF

205 FORMAT('         <kPointList name="',a,'" count="',i0,'">')
    WRITE(fh,205) adjustl(trim(kpts%name)),kpts%nkpt
    IF (kpts%numSpecialPoints<2) THEN
       DO n=1,kpts%nkpt
206       FORMAT('            <kPoint weight="',f12.6,'">',f12.6,' ',f12.6,' ',f12.6,'</kPoint>') 
          WRITE (fh,206) kpts%wtkpt(n), kpts%bk(:,n)
       END DO
       IF (kpts%ntet>0) THEN
          WRITE(fh,207) kpts%ntet
207       FORMAT('            <tetraeder ntet="',i0,'">')
          DO n=1,kpts%ntet
208          FORMAT('          <tet vol="',f12.6,'">',i0,' ',i0,' ',i0,' ',i0,'</tet>')
             WRITE(fh,208) kpts%voltet(n),kpts%ntetra(:,n)
          END DO
          WRITE(fh,'(a)') '            </tetraeder>'
       ELSE
          DO n = 1, kpts%numSpecialPoints
             WRITE(fh,209) TRIM(ADJUSTL(kpts%specialPointNames(n))),kpts%specialPoints(:,n)
209          FORMAT('            <specialPoint name="',a,'">', f10.6,' ',f10.6,' ',f10.6,'</specialPoint>')
          END DO
       END IF
    END IF
    WRITE (fh,'(a)')('         </kPointList>')
    IF (PRESENT(filename)) CLOSE(fh)
  END SUBROUTINE print_xml


  SUBROUTINE add_special_line(kpts,point,name)
    CLASS(t_kpts),INTENT(inout):: kpts
    CHARACTER(len=*),INTENT(in):: name
    REAL,INTENT(in)            :: point(3)

    CHARACTER(len=50),ALLOCATABLE:: names(:)
    REAL,ALLOCATABLE             :: points(:,:)

    IF (kpts%numspecialpoints>0) THEN
       CALL MOVE_ALLOC(kpts%specialPointNames,names)
       CALL MOVE_ALLOC(kpts%specialPoints,points)
       DEALLOCATE(kpts%specialpointindices)
       ALLOCATE(kpts%specialPoints(3,SIZE(names)+1))
       ALLOCATE(kpts%specialPointIndices(SIZE(names)+1))
       ALLOCATE(kpts%specialPointNames(SIZE(names)+1))
       kpts%specialPoints(:,:SIZE(names))=points
       kpts%specialPointNames(:SIZE(names))=names
    ELSE
       ALLOCATE(kpts%specialPoints(3,1))
       ALLOCATE(kpts%specialPointIndices(1))
       ALLOCATE(kpts%specialPointNames(1))
    ENDIF
    kpts%numspecialpoints=kpts%numspecialpoints+1
    kpts%specialPoints(:,kpts%numspecialpoints)=point
    kpts%specialPointNames(kpts%numspecialpoints)=name
  END SUBROUTINE add_special_line

  SUBROUTINE add_special_points_default(kpts,film,cell)
    USE m_judft
    USE m_bravais
    USE m_types_cell
    TYPE(t_kpts),INTENT(inout):: kpts
    LOGICAL,INTENT(in)        :: film
    TYPE(t_cell),INTENT(in)   :: cell

    REAL, PARAMETER :: f12 = 1./2., f14 = 1./4., zro = 0.0
    REAL, PARAMETER :: f34 = 3./4., f38 = 3./8., one = 1.0
    REAL, PARAMETER :: f13 = 1./3., f23 = 2./3.

    INTEGER:: idsyst,idtype
    CALL bravais(cell%amat,idsyst,idtype) 

    IF (.NOT.film) THEN
       IF ( (idsyst == 1).AND.(idtype ==  3) ) THEN       ! fcc
          CALL kpts%add_special_line((/f12,f12,one/) ,"X")
          CALL kpts%add_special_line((/f38,f38,f34/) ,"K")
          CALL kpts%add_special_line((/zro,zro,zro/) ,"g")
          CALL kpts%add_special_line((/f12,f12,f12/) ,"L")
          CALL kpts%add_special_line((/f12,f14,f34/) ,"W")
          CALL kpts%add_special_line((/f12,zro,f12/) ,"X")
          CALL kpts%add_special_line((/zro,zro,zro/) ,"g")
       ENDIF
       IF ( (idsyst == 5).AND.(idtype ==  1) ) THEN       ! rhombohedric (trigonal)
          CALL kpts%add_special_line((/f12,f12, f12/) ,"Z")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/f14,f14,-f14/) ,"K")
          CALL kpts%add_special_line((/f12,f12,-f12/) ,"Z")
          CALL kpts%add_special_line((/f14,f12,-f14/) ,"W")
          CALL kpts%add_special_line((/zro,f12, zro/) ,"L")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/f12,f12, zro/) ,"F")
       ENDIF
       IF ( (idsyst == 4).AND.(idtype ==  1) ) THEN       ! hexagonal
          IF (cell%bmat(1,1)*cell%bmat(2,1)+cell%bmat(1,2)*cell%bmat(2,2) > 0.0) THEN
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/zro,f12, zro/) ,"M")
             CALL kpts%add_special_line((/f13,f13, zro/) ,"K")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
             CALL kpts%add_special_line((/zro,f12, f12/) ,"L")
             CALL kpts%add_special_line((/f13,f13, f12/) ,"H")
             CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
          ELSE                                             ! hexagonal (angle = 60)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
             CALL kpts%add_special_line((/f13,f23, zro/) ,"K")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
             CALL kpts%add_special_line((/f12,f12, f12/) ,"L")
             CALL kpts%add_special_line((/f13,f23, f12/) ,"H")
             CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
          ENDIF
       ENDIF
       IF ( (idsyst == 1).AND.(idtype ==  1) ) THEN       ! simple cubic
          CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/f12,zro, zro/) ,"X")
          CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
          CALL kpts%add_special_line((/f12,f12, f12/) ,"R")
          CALL kpts%add_special_line((/f12,zro, zro/) ,"X")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/f12,f12, f12/) ,"R")
       ENDIF
       IF ( (idsyst == 1).AND.(idtype ==  2) ) THEN       ! body centered cubic
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/f12,-f12,f12/) ,"H")
          CALL kpts%add_special_line((/zro,zro, f12/) ,"N")
          CALL kpts%add_special_line((/f14,f14, f14/) ,"P")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/zro,zro, f12/) ,"N")
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a > c)
          CALL kpts%add_special_line((/f12,f12,-f12/) ,"Z")    ! via Lambda and V)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
          CALL kpts%add_special_line((/-f12,f12,f12/) ,"Z")    ! via Y)
          CALL kpts%add_special_line((/zro,zro, f12/) ,"X")    ! via Delta)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    
          CALL kpts%add_special_line((/zro,f12, zro/) ,"N")    ! via Q)
          CALL kpts%add_special_line((/f14,f14, f14/) ,"P")    ! via W)
          CALL kpts%add_special_line((/zro,zro, f12/) ,"X")
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a < c)
          CALL kpts%add_special_line((/-f12,f12,f12/) ,"Z")    ! via F and Sigma)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
          CALL kpts%add_special_line((/zro,zro, f12/) ,"X")    ! via W)
          CALL kpts%add_special_line((/f14,f14, f14/) ,"P")    ! via Q)
          CALL kpts%add_special_line((/zro,f12, zro/) ,"N")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
          CALL kpts%add_special_line((/f12,f12,-f12/) ,"Z")    ! via U and Y)
          CALL kpts%add_special_line((/f12,f12, zro/) ,"X")
          CALL kpts%add_special_line((/f14,f14, f14/) ,"P")
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
          CALL kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via Y)
          CALL kpts%add_special_line((/f12,f12, zro/) ,"M")    ! via Sigma)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
          CALL kpts%add_special_line((/zro,zro, f12/) ,"Z")    ! via U)
          CALL kpts%add_special_line((/f12,zro, f12/) ,"R")    ! via T)
          CALL kpts%add_special_line((/f12,f12, f12/) ,"A")    ! via S)
          CALL kpts%add_special_line((/zro,zro, f12/) ,"Z")
       ENDIF
       IF ( (idsyst == 3).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
          CALL kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via D)
          CALL kpts%add_special_line((/f12,f12, zro/) ,"S")    ! via C)
          CALL kpts%add_special_line((/zro,f12, zro/) ,"Y")    ! via Delta)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
          CALL kpts%add_special_line((/zro,zro, f12/) ,"Z")    ! via A)
          CALL kpts%add_special_line((/f12,zro, f12/) ,"U")    ! via P)
          CALL kpts%add_special_line((/f12,f12, f12/) ,"R")    ! via E)
          CALL kpts%add_special_line((/zro,f12, f12/) ,"T")    ! via B)
          CALL kpts%add_special_line((/zro,zro, f12/), "Z")
       ENDIF
    ELSE
       WRITE(*,*) 'Note:'
       WRITE(*,*) 'Default k point paths for film band structures'
       WRITE(*,*) 'are experimental. If the generated k point path'
       WRITE(*,*) 'is not correct please specify it directly.'
       IF ( (idsyst == 5).AND.(idtype ==  1) ) THEN       ! rhombohedric (trigonal)
          CALL kpts%add_special_line((/zro,f12, zro/) ,"L")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/f12,f12, zro/) ,"F")
       ENDIF
       IF ( (idsyst == 4).AND.(idtype ==  1) ) THEN       ! hexagonal
          IF (cell%bmat(1,1)*cell%bmat(2,1)+cell%bmat(1,2)*cell%bmat(2,2) > 0.0) THEN
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/zro,f12, zro/) ,"M")
             CALL kpts%add_special_line((/f13,f13, zro/) ,"K")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")

          ELSE                                             ! hexagonal (angle = 60)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
             CALL kpts%add_special_line((/f13,f23, zro/) ,"K")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          ENDIF
       ENDIF
       IF ( (idsyst == 1).AND.(idtype ==  1) ) THEN       ! simple cubic
          CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
          CALL kpts%add_special_line((/f12,zro, zro/) ,"X")
          CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
          CALL kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via Y)
          CALL kpts%add_special_line((/f12,f12, zro/) ,"M")    ! via Sigma)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
       ENDIF
       IF ( (idsyst == 3).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
          CALL kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via D)
          CALL kpts%add_special_line((/f12,f12, zro/) ,"S")    ! via C)
          CALL kpts%add_special_line((/zro,f12, zro/) ,"Y")    ! via Delta)
          CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
       ENDIF
    END IF
    IF (kpts%numspecialPoints<2) CALL judft_error("Not enough special points given and no default found")
  END SUBROUTINE add_special_points_default




END MODULE m_types_kpts
