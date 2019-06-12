!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts
  PRIVATE
  TYPE t_kpts
     INTEGER               :: nkpt=0
     INTEGER               :: ntet=0
     LOGICAL               :: l_gamma=.false.
     !(3,nkpt) k-vectors internal units
     REAL,ALLOCATABLE      :: bk(:,:)
     !(nkpts) weights
     REAL,ALLOCATABLE      :: wtkpt(:)
     INTEGER               :: nkptf=0   !<k-vectors in full BZ
     INTEGER               :: nkpt3(3)=[0,0,0]
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
   contains
     procedure :: init_defaults
     procedure :: init_by_density
     procedure :: init_by_number
     procedure :: init_by_grid
     procedure :: init_special
     !procedure :: read_xml
     procedure :: add_special_line
  ENDTYPE t_kpts

  public :: t_kpts
CONTAINS

  
  SUBROUTINE add_special_line(kpts,point,name)
    CLASS(t_kpts),INTENT(inout):: kpts
    CHARACTER(len=*),INTENT(in):: name
    REAL,INTENT(in)            :: point

    IF (kpts%numspecialpoints>0) THEN
       CALL move_ALLOC(kpts%specialPointNames,names)
       CALL move_ALLOC(kpts%specialPoints,points)
       DEALLOCATE(kpts%specialpointindices)
       ALLOCATE(kpts%specialPoints(3,SIZE(names)+1))
       ALLOCATE(kpts%specialPointIndices(SIZE(names)+1))
       ALLOCATE(kpts%specialPointNames(SIZE(names)+1))
       kpts%specialPoints(3,:SIZE(names))=points
       kpts%specialPointNames(:SIZE(names))=names
    ELSE
       ALLOCATE(kpts%specialPoints(3,1))
       ALLOCATE(kpts%specialPointIndices(1))
       ALLOCATE(kpts%specialPointNames(1))
    ENDIF
    kpts%numspecialpoints=kpts%numspecialpoints+1
    kpts%specialPoints(3,kpts%numspecialpoints)=point
    kpts%specialPointNames(kpts%numspecialpoints)=name
  END SUBROUTINE add_special_line
  
  SUBROUTINE init_special(kpts,cell)
    CLASS(t_kpts),INTENT(inout):: kpts
    LOGICAL,INTENT(IN)         :: film
    TYPE(t_cell),intent(IN)    :: cell

    IF (kpts%numSpecialPoints<2) CALL get_special_points(film,cell)
    kpts%nkpt=MAX(kpts%nkpt,kpts%numSpecialPoints)
    !all sepecial kpoints are now set already

    !Distances
    lastp=0
    DO i=1,kpts%numSpecialPoints
       nextp(:)=MATMUL(kpts%specialPoints(:,i),cell%bmat)
       d(i)=SQRT(dot_PRODUCT(nextp-lastp,nextp-lastp))
       lastp=nextp
    ENDDO
    d(1)=0.0
    !Distribute points
    nk(1)=0
    DO i=2,kpts%numSpecialPoints
       nk(i)=NINT((kpts%nkpt-kpts%numSpecialPoints)*(d(i)/SUM(d)))
    ENDDO

    ALLOCATE(kpts%bk(kpts%numSpecialPoints+SUM(nk)))
    
    !Generate lines
    kpts%nkpt=1
    DO i=1,kpts%numSpecialPoints-1
       kpts%bk(:,kpts%nkpt)=kpts%specialPoints(:,i)
       kpts%specialPointIndices(i)=kpts%nkpt
       kpts%nkpt=kpts%nkpt+1
       d=(kpts%specialPoints(:,i+1)-kpts%specialPoints(:,i))/(nk(i)+2) 
       DO ii=1,nk(i)
          kpts%bk(:,kpts%nkpt)=kpts%specialPoints(:,i)+d*ii
          kpts%nkpt=kpts%nkpt+1
       ENDDO
    ENDDO
    kpts%bk(:,kpts%nkpt)=kpts%specialPoints(:,kpts%numSpecialPoints)
    kpts%specialPointIndices(kpts%numSpecialPoints)=kpts%nkpt
    ALLOCATE(kpts%wtkpt(kpts%nkpt))
    kpts%wtkpt = 1.0
  END SUBROUTINE init_special
    
          
  SUBROUTINE init_defaults(kpts,film,area,cell,sym)
    class(t_kpts),intent(out):: kpts
    LOGICAL,INTENT(in)       :: film
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym

    integer:: nkpt
    IF (film) THEN
       nkpt = MAX(nint((3600/cell%area)/sym%nop2),1)
    ELSE
       nkpt = MAX(nint((216000/cell%omtil)/sym%nop),1)
    ENDIF
    call kpts%init_by_number(nkpt)
  end subroutine init_defaults

  SUBROUTINE init_by_density(kpts,density,cell,sym)
    CLASS(t_kpts),INTENT(out):: kpts
    REAL,INTENT(in)          :: density
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym
    REAL    :: length
    INTEGER :: n,grid(3)

    do n=1,3
       length=sqrt(dot_product(cell%bmat(n,:),cell%bmat(n,:)))  !TODO why not bmat(:,n)???
       grid(n)=ceiling(density*length)
    end do
    call kpts%init_by_grid(grid)
    
  end subroutine init_by_density

  SUBROUTINE init_by_number(kpts,nkpt,film,cell,sym)
    USE m_divi
    CLASS(t_kpts),INTENT(out):: kpts
    INTEGER,INTENT(IN)       :: nkpt
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym
    
    INTEGER :: grid(3)

    CALL divi(nkpt,cell%bmat,film,sym%nop,syom%nop2,grid)
    CALL kpts%init_by_grid(grid)
  END SUBROUTINE init_by_number

  SUBROUTINE init_by_grid(kpts,sym,cell,grid,film,tria,l_soc_or_ss)
    !----------------------------------------------------------------------+
    ! Generate a k-point file with approx. nkpt k-pts or a Monkhorst-Pack  |
    ! set with nmod(i) divisions in i=x,y,z direction. Interface to kptmop |
    ! and kpttet routines of the MD-programm.                              |
    !                                                          G.B. 07/01  |
    !----------------------------------------------------------------------+
    USE m_constants
    USE m_bravais
    USE m_brzone2
    USE m_kptmop
    USE m_kpttet
    USE m_gen_bz
    USE m_types

    IMPLICIT NONE
    CLASS(t_kpts),INTENT(out):: kpts

    TYPE(t_sym),     INTENT(IN)    :: sym
    TYPE(t_cell),    INTENT(IN)    :: cell
    INTEGER,INTENT(IN)             :: grid(3)
    LOGICAL,INTENT(IN)             :: film,tria,l_soc_or_ss


    INTEGER, PARAMETER :: nop48  = 48
    INTEGER, PARAMETER :: mface  = 51
    INTEGER, PARAMETER :: mdir   = 10
    INTEGER, PARAMETER :: nbsz   =  3
    INTEGER, PARAMETER :: ibfile =  6
    INTEGER, PARAMETER :: nv48   = (2*nbsz+1)**3+48

    INTEGER ndiv3              ! max. number of tetrahedrons (< 6*(kpts%nkpt+1)

    REAL, ALLOCATABLE    :: vkxyz(:,:)  ! vector of kpoint generated; in cartesian representation
    REAL, ALLOCATABLE    :: wghtkp(:)   !   associated with k-points for BZ integration
    INTEGER, ALLOCATABLE :: ntetra(:,:) ! corners of the tetrahedrons
    REAL, ALLOCATABLE    :: voltet(:)   ! voulmes of the tetrahedrons
    REAL, ALLOCATABLE    :: vktet(:,:)

    REAL    divis(4)           ! Used to find more accurate representation of k-points
    ! vklmn(i,kpt)/divis(i) and weights as wght(kpt)/divis(4)
    REAL    bltv(3,3)          ! cartesian Bravais lattice basis (a.u.)
    REAL    rltv(3,3)          ! reciprocal lattice basis (2\pi/a.u.)
    REAL    ccr(3,3,nop48)     ! rotation matrices in cartesian repr.
    REAL    rlsymr(3,3,nop48)  ! rotation matrices in reciprocal lattice basis representation
    REAL    talfa(3,nop48)     ! translation vector associated with (non-symmorphic)
    ! symmetry elements in Bravais lattice representation
    INTEGER ncorn,nedge,nface  ! number of corners, faces and edges of the IBZ
    REAL    fnorm(3,mface)     ! normal vector of the planes bordering the IBZ
    REAL    fdist(mface)       ! distance vector of the planes bordering t IBZ
    REAL    cpoint(3,mface)    ! cartesian coordinates of corner points of IBZ
    REAL    xvec(3)            ! arbitrary vector lying in the IBZ

    INTEGER idsyst   ! crystal system identification in MDDFT programs
    INTEGER idtype   ! lattice type identification in MDDFT programs

    INTEGER idimens  ! number of dimensions for k-point set (2 or 3)
    INTEGER nreg     ! 1 kpoints in full BZ; 0 kpoints in irrBZ
    INTEGER nfulst   ! 1 kpoints ordered in full stars
    !    (meaningful only for nreg =1; full BZ)
    INTEGER nbound   ! 0 no primary points on BZ boundary;
    ! 1 with boundary points (not for BZ integration!!!)
    INTEGER ikzero   ! 0 no shift of k-points;
    ! 1 shift of k-points for better use of sym in irrBZ
    REAL    kzero(3) ! shifting vector to bring one k-point to or 
    ! away from (0,0,0) (for even/odd nkpt3)

    INTEGER i,j,k,l,mkpt,addSym,nsym
    LOGICAL random,trias
    REAL help(3),binv(3,3),rlsymr1(3,3),ccr1(3,3)

    IF (ANY(grid==0)) THEN
       PRINT *,"Warning, k-point grid dimension increased to 1"
       WHERE(grid==0) grid=1
    END IF


    IF (l_gamma) THEN
       IF (tria) CALL judft_error("tria and l_gamma incompatible")
       CALL kptgen_hybrid(film,grid,cell,sym,kpts,l_soc_or_ss)
    ELSE
       !------------------------------------------------------------
       !
       !        idsyst         idtype 
       !
       !   1  cubic          primitive
       !   2  tetragonal     body centered
       !   3  orthorhombic   face centered
       !   4  hexagonal      A-face centered
       !   5  trigonal       B-face centered
       !   6  monoclinic     C-face centered
       !   7  triclinic 
       !
       ! --->   for 2 dimensions only the following Bravais lattices exist:
       !
       !    TYPE                    EQUIVALENT 3-DIM        idsyst/idtype
       !   square               = p-tetragonal ( 1+2 axis )      2/1
       !   rectangular          = p-orthorhomb ( 1+2 axis )      3/1
       !   centered rectangular = c-face-orthorhomb( 1+2 axis)   3/6
       !   hexagonal            = p-hexagonal  ( 1+2 axis )      4/1
       !   oblique              = p-monoclinic ( 1+2 axis )      6/1
       !
       !------------------------------------------------------------


       CALL bravais(cell%amat,idsyst,idtype) 

       nsym = MERGE(sym%nop2,sym%nop,film)    
       nbound  = MERGE(1,0,film.AND.tria)
       random  = tria.AND..NOT.film
       idimens = MERGE(2,3,film)

       ! Lattice information

       bltv=TRANSPOSE(cell%amat)
       binv=TRANSPOSE(cell%bmat)/tpi_conts
       rltv=TRANSPOSE(cell%bmat)
       DO i=1,nsym
          rlsymr(:,:,i)=REAL(TRANSPOSE(sym%mat(:,:,i))
       ENDDO

       talfa=MATMUL(bltv,sym%tau(:,:nsym))
       DO i = 1, nsym
          ccr(:,:,i) = MATMUL(MATMUL(binv(:,:),rlsymr(:,:,i)),bltv(:,:))
       END DO
       DO i = 1, nsym
          rlsymr(:,:,i)=TRANSPOSE(rlsymr(:,:,i))
          ccr(:,:,i)=TRANSPOSE(ccr(:,:,i))
       END DO

       IF ((.not.l_soc_or_ss).AND.(2*nsym<nop48)) THEN
          IF ((film.AND.(.not.sym%invs2)).OR.((.not.film).AND.(.not.sym%invs))) THEN
             addSym = 0
             ! Note: We have to add the negative of each symmetry operation
             !       to exploit time reversal symmetry. However, if the new
             !       symmetry operation is the identity matrix it is excluded.
             !       This is the case iff it is (-Id) + a translation vector.
             DO i = 1, nsym
                ! This test assumes that ccr(:,:,1) is the identity matrix.
                IF(.NOT.ALL(ABS(ccr(:,:,1)+ccr(:,:,i)).LT.10e-10) ) THEN
                   ccr(:,:,nsym+addSym+1 ) = -ccr(:,:,i)
                   rlsymr(:,:,nsym+addSym+1 ) = -rlsymr(:,:,i)
                   addSym = addSym + 1
                END IF
             END DO
             nsym = nsym + addSym
          END IF
       END IF

       ! brzone and brzone2 find the corner-points, the edges, and the
       ! faces of the irreducible wedge of the brillouin zone (IBZ).
       CALL brzone2(rltv,nsym,ccr,mface,nbsz,nv48,cpoint,xvec,ncorn,nedge,nface,fnorm,fdist)

       IF (nbound.EQ.1) THEN
          mkpt = PRODUCT(2*grid(:idimens)+1))
       ELSE
          mkpt=PRODUCT(grid(:idimens))
       END IF
       ALLOCATE (vkxyz(3,mkpt),wghtkp(mkpt) )


       IF (tria.AND.random) THEN
          ! Calculate the points for tetrahedron method
          ndiv3 = 6*(mkpt+1)
          ALLOCATE (voltet(ndiv3),vktet(3,mkpt),ntetra(4,ndiv3))
          kpts%nkpt=mkpt
          CALL kpttet(0,mkpt,ndiv3,&
               rltv,cell%omtil,nsym,ccr,mdir,mface,&
               ncorn,nface,fdist,fnorm,cpoint,voltet,ntetra,kpts%ntet,vktet,&
               kpts%nkpt,vkxyz,wghtkp)
       ELSE
          ! Now calculate Monkhorst-Pack k-points:
          CALL kptmop(idsyst,idtype,grid,&
               rltv,bltv,nbound,idimens,xvec,fnorm,fdist,ncorn,nface,&
               nedge,cpoint,nsym,ccr,rlsymr,talfa,mkpt,mface,mdir,&
               kpts%nkpt,vkxyz,wghtkp)
       END IF

       DO j=1,kpts%nkpt
          vkxyz(:,j)=MATMUL(vkxyz(:,j),cell%amat)/tpi_const
       END DO

       ALLOCATE(kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt))
       kpts%bk(:,:) = vkxyz(:,:kpts%nkpt)
       kpts%wtkpt(:) = wghtkp(:kpts%nkpt)

       IF (tria.AND.random) THEN
          ALLOCATE(kpts%ntetra(4,kpts%ntet))
          ALLOCATE(kpts%voltet(kpts%ntet))
          DO j = 1, kpts%ntet
             kpts%ntetra(:,j) = ntetra(:,j)
             kpts%voltet(j) = ABS(voltet(j))
          END DO
       END IF
    ENDIF
    !Generate kpts in full BZ
    IF (.NOT.tria) CALL gen_bz(kpts,sym)


  END SUBROUTINE init_by_grid


SUBROUTINE add_special_points_default(kpts,film,amat)
    TYPE(t_kpts),INTENT(inout):: kpts
    LOGICAL,INTENT(in)        :: film
    REAL,INTENT(in)           :: amat(3,3)
    
    REAL, PARAMETER :: f12 = 1./2., f14 = 1./4., zro = 0.0
    REAL, PARAMETER :: f34 = 3./4., f38 = 3./8., one = 1.0
    REAL, PARAMETER :: f13 = 1./3., f23 = 2./3.

    INTEGER:: idsyst,idtype
    CALL bravais(amat,idsyst,idtype) 
    
      IF (.NOT.film) THEN
         IF ( (idsyst == 1).AND.(idtype ==  3) ) THEN       ! fcc
            CALL kpts%add_special_line((/f12,f12,one/) ,"X")
            call kpts%add_special_line((/f38,f38,f34/) ,"K")
            call kpts%add_special_line((/zro,zro,zro/) ,"g")
            call kpts%add_special_line((/f12,f12,f12/) ,"L")
            call kpts%add_special_line((/f12,f14,f34/) ,"W")
            call kpts%add_special_line((/f12,zro,f12/) ,"X")
            call kpts%add_special_line((/zro,zro,zro/) ,"g")
         ENDIF
         IF ( (idsyst == 5).AND.(idtype ==  1) ) THEN       ! rhombohedric (trigonal)
            CALL kpts%add_special_line((/f12,f12, f12/) ,"Z")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/f14,f14,-f14/) ,"K")
            call kpts%add_special_line((/f12,f12,-f12/) ,"Z")
            call kpts%add_special_line((/f14,f12,-f14/) ,"W")
            call kpts%add_special_line((/zro,f12, zro/) ,"L")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/f12,f12, zro/) ,"F")
         ENDIF
         IF ( (idsyst == 4).AND.(idtype ==  1) ) THEN       ! hexagonal
            IF (bmat(1,1)*bmat(2,1)+bmat(1,2)*bmat(2,2) > 0.0) THEN
               call kpts%add_special_line((/zro,zro, zro/) ,"g")
               call kpts%add_special_line((/zro,f12, zro/) ,"M")
               call kpts%add_special_line((/f13,f13, zro/) ,"K")
               call kpts%add_special_line((/zro,zro, zro/) ,"g")
               call kpts%add_special_line((/zro,zro, f12/) ,"A")
               call kpts%add_special_line((/zro,f12, f12/) ,"L")
               call kpts%add_special_line((/f13,f13, f12/) ,"H")
               call kpts%add_special_line((/zro,zro, f12/) ,"A")
            ELSE                                             ! hexagonal (angle = 60)
               call kpts%add_special_line((/zro,zro, zro/) ,"g")
               call kpts%add_special_line((/f12,f12, zro/) ,"M")
               call kpts%add_special_line((/f13,f23, zro/) ,"K")
               call kpts%add_special_line((/zro,zro, zro/) ,"g")
               call kpts%add_special_line((/zro,zro, f12/) ,"A")
               call kpts%add_special_line((/f12,f12, f12/) ,"L")
               call kpts%add_special_line((/f13,f23, f12/) ,"H")
               call kpts%add_special_line((/zro,zro, f12/) ,"A")
            ENDIF
         ENDIF
         IF ( (idsyst == 1).AND.(idtype ==  1) ) THEN       ! simple cubic
            CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/f12,zro, zro/) ,"X")
            call kpts%add_special_line((/f12,f12, zro/) ,"M")
            call kpts%add_special_line((/f12,f12, f12/) ,"R")
            call kpts%add_special_line((/f12,zro, zro/) ,"X")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/f12,f12, f12/) ,"R")
         ENDIF
         IF ( (idsyst == 1).AND.(idtype ==  2) ) THEN       ! body centered cubic
            CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/f12,-f12,f12/) ,"H")
            call kpts%add_special_line((/zro,zro, f12/) ,"N")
            call kpts%add_special_line((/f14,f14, f14/) ,"P")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/zro,zro, f12/) ,"N")
         ENDIF
         IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a > c)
            call kpts%add_special_line((/f12,f12,-f12/) ,"Z")    ! via Lambda and V)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
            call kpts%add_special_line((/-f12,f12,f12/) ,"Z")    ! via Y)
            call kpts%add_special_line((/zro,zro, f12/) ,"X")    ! via Delta)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    
            call kpts%add_special_line((/zro,f12, zro/) ,"N")    ! via Q)
            call kpts%add_special_line((/f14,f14, f14/) ,"P")    ! via W)
            call kpts%add_special_line((/zro,zro, f12/) ,"X")
         ENDIF
         IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a < c)
            CALL kpts%add_special_line((/-f12,f12,f12/) ,"Z")    ! via F and Sigma)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
            call kpts%add_special_line((/zro,zro, f12/) ,"X")    ! via W)
            call kpts%add_special_line((/f14,f14, f14/) ,"P")    ! via Q)
            call kpts%add_special_line((/zro,f12, zro/) ,"N")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
            call kpts%add_special_line((/f12,f12,-f12/) ,"Z")    ! via U and Y)
            call kpts%add_special_line((/f12,f12, zro/) ,"X")
            call kpts%add_special_line((/f14,f14, f14/) ,"P")
         ENDIF
         IF ( (idsyst == 2).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
            call kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via Y)
            call kpts%add_special_line((/f12,f12, zro/) ,"M")    ! via Sigma)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
            call kpts%add_special_line((/zro,zro, f12/) ,"Z")    ! via U)
            call kpts%add_special_line((/f12,zro, f12/) ,"R")    ! via T)
            call kpts%add_special_line((/f12,f12, f12/) ,"A")    ! via S)
            call kpts%add_special_line((/zro,zro, f12/) ,"Z")
         ENDIF
         IF ( (idsyst == 3).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
            call kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via D)
            call kpts%add_special_line((/f12,f12, zro/) ,"S")    ! via C)
            call kpts%add_special_line((/zro,f12, zro/) ,"Y")    ! via Delta)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
            call kpts%add_special_line((/zro,zro, f12/) ,"Z")    ! via A)
            call kpts%add_special_line((/f12,zro, f12/) ,"U")    ! via P)
            call kpts%add_special_line((/f12,f12, f12/) ,"R")    ! via E)
            call kpts%add_special_line((/zro,f12, f12/) ,"T")    ! via B)
            syp(:,10)= (/zro,zro, f12/)  ; ssy(8) = "Z"
         ENDIF
      ELSE
         WRITE(*,*) 'Note:'
         WRITE(*,*) 'Default k point paths for film band structures'
         WRITE(*,*) 'are experimental. If the generated k point path'
         WRITE(*,*) 'is not correct please specify it directly.'
         IF ( (idsyst == 5).AND.(idtype ==  1) ) THEN       ! rhombohedric (trigonal)
            call kpts%add_special_line((/zro,f12, zro/) ,"L")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/f12,f12, zro/) ,"F")
         ENDIF
         IF ( (idsyst == 4).AND.(idtype ==  1) ) THEN       ! hexagonal
            IF (bmat(1,1)*bmat(2,1)+bmat(1,2)*bmat(2,2) > 0.0) THEN
               call kpts%add_special_line((/zro,zro, zro/) ,"g")
               call kpts%add_special_line((/zro,f12, zro/) ,"M")
               call kpts%add_special_line((/f13,f13, zro/) ,"K")
               call kpts%add_special_line((/zro,zro, zro/) ,"g")

            ELSE                                             ! hexagonal (angle = 60)
               call kpts%add_special_line((/zro,zro, zro/) ,"g")
               call kpts%add_special_line((/f12,f12, zro/) ,"M")
               call kpts%add_special_line((/f13,f23, zro/) ,"K")
               call kpts%add_special_line((/zro,zro, zro/) ,"g")
            ENDIF
         ENDIF
         IF ( (idsyst == 1).AND.(idtype ==  1) ) THEN       ! simple cubic
            call kpts%add_special_line((/f12,f12, zro/) ,"M")
            call kpts%add_special_line((/zro,zro, zro/) ,"g")
            call kpts%add_special_line((/f12,zro, zro/) ,"X")
            call kpts%add_special_line((/f12,f12, zro/) ,"M")
         ENDIF
         IF ( (idsyst == 2).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
            call kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via Y)
            call kpts%add_special_line((/f12,f12, zro/) ,"M")    ! via Sigma)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
         ENDIF
         IF ( (idsyst == 3).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
            call kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via D)
            call kpts%add_special_line((/f12,f12, zro/) ,"S")    ! via C)
            call kpts%add_special_line((/zro,f12, zro/) ,"Y")    ! via Delta)
            call kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
         ENDIF
      END IF
      IF (kpts%numspecialPoints<2) CALL judft_error("Not enough special points given and no default found")
    END SUBROUTINE add_special_points_default
  

  
END MODULE m_types_kpts
