!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_make_kpoints
  USE m_juDFT
  use m_types_kpts
  IMPLICIT NONE
  private
  public :: make_kpoints
CONTAINS
  SUBROUTINE make_kpoints(kpts,cell,sym,film,l_socorss,str)
    USE m_types_kpts
    USE m_types_cell
    USE m_types_sym
    TYPE(t_kpts),INTENT(out)::kpts
    TYPE(t_cell),INTENT(in) ::cell
    TYPE(t_sym),INTENT(in)  ::sym
    LOGICAL,INTENT(in)::l_socorss,film
    CHARACTER(len=*),INTENT(inout)::str

    LOGICAL:: tria,l_gamma,l_soc_or_ss
    REAL   :: den
    INTEGER:: nk,grid(3)
    character(len=20)::name=""
    !defaults
    l_soc_or_ss=l_socorss
    tria=.false.;l_gamma=.false.

    IF (judft_was_argument("-k")) THEN
       IF (LEN_TRIM(str)>1) CALL judft_error("Do not specify k-points in file and on command line")
       str=judft_string_for_argument("-k")
    END IF

    !set name
    IF (INDEX(str,"#")>0) THEN
       name=str(:INDEX(str,"#")-1)
       str=trim(adjustl(str(INDEX(str,"#")+1:)))
    END IF
    str=ADJUSTL(str)
    DO WHILE(INDEX(str,'@')>0)
       IF (INDEX(str,'tria@')==1) THEN
          tria=.TRUE.
          str=str(6:)
       ENDIF
       IF (INDEX(str,'gamma@')==1) THEN
          l_gamma=.TRUE.
          str=str(7:)
       ENDIF
       IF (INDEX(str,'soc@')==1) THEN
          l_soc_or_ss=.TRUE.
          str=str(5:)
       ENDIF
    END DO
    IF (INDEX(str,'den=')==1) THEN
       str=str(5:)
       READ(str,*) den
       CALL init_by_density(kpts,den,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSEIF(INDEX(str,'nk=')==1) THEN
       str=str(4:)
       READ(str,*) nk
       CALL init_by_number(kpts,nk,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSEIF(INDEX(str,'grid=')==1) THEN
       str=str(6:)
       READ(str,*) grid
       CALL init_by_grid(kpts,grid,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSEIF(INDEX(str,'file')==1) THEN
       CALL init_by_kptsfile(kpts,film)
    ELSEIF(LEN_TRIM(str)<1.OR.INDEX(ADJUSTL(str),'#')==1) THEN
       CALL init_defaults(kpts,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSE
       CALL judft_error(("Could not process -k argument:"//str))
    ENDIF

    if (len_trim(name)>0) kpts%name=name
  END SUBROUTINE make_kpoints


  SUBROUTINE init(kpts,cell,sym,film)
    USE m_types_cell
    USE m_types_sym
    CLASS(t_kpts),INTENT(inout):: kpts
    TYPE(t_cell),INTENT(IN)    :: cell
    TYPE(t_sym),INTENT(IN)     :: sym
    LOGICAL,INTENT(IN)         :: film

    INTEGER :: n

    IF (kpts%numSpecialPoints>2) THEN
       CALL init_special(kpts,cell,film)
       return
    ENDIF
    DO n=1,kpts%nkpt
       kpts%l_gamma=kpts%l_gamma.OR.ALL(ABS(kpts%bk(:,n))<1E-9)
    ENDDO
    IF (kpts%nkptf==0) CALL gen_bz(kpts,sym)
  END SUBROUTINE init

  SUBROUTINE init_by_kptsfile(kpts,film)
    CLASS(t_kpts),INTENT(out):: kpts
    LOGICAL,INTENT(in)       :: film

    LOGICAL :: l_exist
    REAL    :: scale,wscale
    INTEGER :: n,ios

    INQUIRE(file='kpts',exist=l_exist)
    IF (.NOT.l_exist) CALL judft_error("Could not read 'kpts' file")
    OPEN(99,file='kpts')
    READ(99,"(i5,2f20.10,3x,l1)",iostat=ios) kpts%nkpt,scale,wscale
    ALLOCATE(kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt))
    DO n=1,kpts%nkpt
       IF (.NOT.film) THEN
          READ(99,*,iostat=ios) kpts%bk(:,n),kpts%wtkpt(n)
       ELSE
          READ(99,*,iostat=ios) kpts%bk(1:2,n),kpts%wtkpt(n)
          kpts%bk(3,n)=0.0
       ENDIF
    ENDDO
    !check for tetraeder
    READ(99,*,err=100,end=100) kpts%ntet
    ALLOCATE(kpts%ntetra(4,kpts%ntet),kpts%voltet(kpts%ntet))
    read(99,*) kpts%ntetra
    read(99,*) kpts%voltet
100 CLOSE(99)
    IF (ios.NE.0) CALL judft_error("Error while reading 'kpts' file")
    IF (scale>0.0) kpts%bk=kpts%bk/scale
    IF (wscale>0.0) kpts%wtkpt=kpts%wtkpt/wscale
  END SUBROUTINE init_by_kptsfile



  SUBROUTINE init_special(kpts,cell,film)
    USE m_types_cell
    CLASS(t_kpts),INTENT(inout):: kpts
    LOGICAL,INTENT(IN)         :: film
    TYPE(t_cell),INTENT(IN)    :: cell

    REAL:: nextp(3),lastp(3),d(MAX(kpts%nkpt,kpts%numSpecialPoints))
    INTEGER:: nk(MAX(kpts%nkpt,kpts%numSpecialPoints)),i,ii
    IF (kpts%numSpecialPoints<2) CALL add_special_points_default(kpts,film,cell)
    kpts%nkpt=MAX(kpts%nkpt,kpts%numSpecialPoints)
    !all sepecial kpoints are now set already

    !Distances
    lastp=0
    DO i=1,kpts%numSpecialPoints
       nextp(:)=MATMUL(kpts%specialPoints(:,i),cell%bmat)
       d(i)=SQRT(DOT_PRODUCT(nextp-lastp,nextp-lastp))
       lastp=nextp
    ENDDO
    d(1)=0.0
    !Distribute points
    nk(1)=0
    DO i=2,kpts%numSpecialPoints
       nk(i)=NINT((kpts%nkpt-kpts%numSpecialPoints)*(d(i)/SUM(d)))
    ENDDO

    ALLOCATE(kpts%bk(3,kpts%numSpecialPoints+SUM(nk)))

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


  SUBROUTINE init_defaults(kpts,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    USE m_types_cell
    USE m_types_sym
    CLASS(t_kpts),INTENT(out):: kpts
    LOGICAL,INTENT(in)       :: film,tria,l_soc_or_ss,l_gamma
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym

    INTEGER:: nkpt
    IF (film) THEN
       nkpt = MAX(NINT((3600/cell%area)/sym%nop2),1)
    ELSE
       nkpt = MAX(NINT((216000/cell%omtil)/sym%nop),1)
    ENDIF
    CALL init_by_number(kpts,nkpt,cell,sym,film,tria,l_soc_or_ss,l_gamma)
  END SUBROUTINE init_defaults

  SUBROUTINE init_by_density(kpts,density,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    USE m_types_cell
    USE m_types_sym
    CLASS(t_kpts),INTENT(out):: kpts
    REAL,INTENT(in)          :: density
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym
    LOGICAL,INTENT(IN)             :: film,tria,l_soc_or_ss,l_gamma
    REAL    :: length
    INTEGER :: n,grid(3)

    DO n=1,3
       length=SQRT(DOT_PRODUCT(cell%bmat(n,:),cell%bmat(n,:)))  !TODO why not bmat(:,n)???
       grid(n)=CEILING(density*length)
    END DO
    CALL init_by_grid(kpts,grid,cell,sym,film,tria,l_soc_or_ss,l_gamma)

  END SUBROUTINE init_by_density

  SUBROUTINE init_by_number(kpts,nkpt,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    USE m_divi
    USE m_types_cell
    USE m_types_sym
    CLASS(t_kpts),INTENT(out):: kpts
    INTEGER,INTENT(IN)       :: nkpt
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym
    LOGICAL,INTENT(IN)             :: film,tria,l_soc_or_ss,l_gamma

    INTEGER :: grid(3)

    CALL divi(nkpt,cell%bmat,film,sym%nop,sym%nop2,grid)
    CALL init_by_grid(kpts,grid,cell,sym,film,tria,l_soc_or_ss,l_gamma)
  END SUBROUTINE init_by_number

  SUBROUTINE init_by_grid(kpts,grid,cell,sym,film,tria,l_soc_or_ss,l_gamma)
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
    USE m_types_cell
    USE m_types_sym
    USE m_kptgen_hybrid
    IMPLICIT NONE
    CLASS(t_kpts),INTENT(out):: kpts

    TYPE(t_sym),     INTENT(IN)    :: sym
    TYPE(t_cell),    INTENT(IN)    :: cell
    INTEGER,INTENT(INout)          :: grid(3)
    LOGICAL,INTENT(IN)             :: film,tria,l_soc_or_ss,l_gamma


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
       binv=TRANSPOSE(cell%bmat)/tpi_const
       rltv=TRANSPOSE(cell%bmat)
       DO i=1,nsym
          rlsymr(:,:,i)=REAL(TRANSPOSE(sym%mrot(:,:,i)))
       ENDDO

       talfa(:,:nsym)=MATMUL(bltv,sym%tau(:,:nsym))
       DO i = 1, nsym
          ccr(:,:,i) = MATMUL(MATMUL(binv(:,:),rlsymr(:,:,i)),bltv(:,:))
       END DO
       DO i = 1, nsym
          rlsymr(:,:,i)=TRANSPOSE(rlsymr(:,:,i))
          ccr(:,:,i)=TRANSPOSE(ccr(:,:,i))
       END DO

       IF ((.NOT.l_soc_or_ss).AND.(2*nsym<nop48)) THEN
          IF ((film.AND.(.NOT.sym%invs2)).OR.((.NOT.film).AND.(.NOT.sym%invs))) THEN
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
          mkpt = PRODUCT((2*grid(:idimens)+1))
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


  SUBROUTINE gen_bz( kpts,sym)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! gen_bz generates the (whole) Brillouin zone from the          !
    ! (irreducible) k-points given                                  !
    !                                     M.Betzinger (09/07)       !
    !                        Refactored in 2017 by G.M.             !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     bk     ::    irreducible k-points
    !     nkpt   ::    number of irr. k-points
    !     bkf    ::    all k-points
    !     nkptf  ::    number of all k-points
    !     bkp    ::    k-point parent
    !     bksym  ::    symmetry operation, that connects the parent
    !                  k-point with the current one
    USE m_juDFT
    USE m_types_sym
    TYPE(t_kpts),INTENT(INOUT) :: kpts
    TYPE(t_sym),INTENT(IN)     :: sym
    !  - local scalars -
    INTEGER                 ::  ic,iop,ikpt,ikpt1
    LOGICAL                 ::  l_found

    !  - local arrays -
    INTEGER,ALLOCATABLE     ::  iarr(:)
    REAL                    ::  rrot(3,3,2*sym%nop),rotkpt(3)
    REAL,ALLOCATABLE        ::  rarr1(:,:)

    INTEGER:: nsym,ID_mat(3,3)


    nsym=sym%nop
    IF (.NOT.sym%invs) nsym=2*sym%nop

    ALLOCATE (kpts%bkf(3,nsym*kpts%nkpt))
    ALLOCATE (kpts%bkp(nsym*kpts%nkpt))
    ALLOCATE (kpts%bksym(nsym*kpts%nkpt))

    ! Generate symmetry operations in reciprocal space
    DO iop=1,nsym
       IF( iop .LE. sym%nop ) THEN
          rrot(:,:,iop) = TRANSPOSE( sym%mrot(:,:,sym%invtab(iop)) )
       ELSE
          rrot(:,:,iop) = -rrot(:,:,iop-sym%nop)
       END IF
    END DO


    !Add existing vectors to list of full vectors
    id_mat=0
    ID_mat(1,1)=1;ID_mat(2,2)=1;ID_mat(3,3)=1
    IF (ANY(sym%mrot(:,:,1).NE.ID_mat)) CALL judft_error("Identity must be first symmetry operation",calledby="gen_bz")

    ic=0
    DO iop=1,nsym
       DO ikpt=1,kpts%nkpt
          rotkpt = MATMUL(rrot(:,:,iop), kpts%bk(:,ikpt))
          !transform back into 1st-BZ (Do not use nint to deal properly with inaccuracies)
          do while(any(rotkpt<-0.5))
             where (rotkpt<-0.5) rotkpt=rotkpt+1.0
          enddo
          do while(any(rotkpt>0.5+1E-8))
             where (rotkpt>0.5+1E-8) rotkpt=rotkpt-1.0
          enddo
          DO ikpt1=1,ic
             IF (MAXVAL(ABS(kpts%bkf(:,ikpt1) - rotkpt)).LE.1e-07) EXIT
          END DO

          IF(ikpt1>ic) THEN !new point
             ic = ic + 1
             kpts%bkf(:,ic) = rotkpt
             kpts%bkp(ic) = ikpt
             kpts%bksym(ic) = iop
          END IF
       END DO
    END DO

    kpts%nkptf = ic
    ! Reallocate bkf, bkp, bksym
    ALLOCATE (iarr(kpts%nkptf))
    iarr = kpts%bkp(:kpts%nkptf)
    DEALLOCATE(kpts%bkp)
    ALLOCATE (kpts%bkp(kpts%nkptf))
    kpts%bkp = iarr
    iarr= kpts%bksym(:kpts%nkptf)
    DEALLOCATE (kpts%bksym )
    ALLOCATE (kpts%bksym(kpts%nkptf))
    kpts%bksym = iarr
    DEALLOCATE(iarr)
    ALLOCATE (rarr1(3,kpts%nkptf))
    rarr1 = kpts%bkf(:,:kpts%nkptf)
    DEALLOCATE (kpts%bkf )
    ALLOCATE (kpts%bkf(3,kpts%nkptf))
    kpts%bkf = rarr1
    DEALLOCATE(rarr1)

  END SUBROUTINE gen_bz



END MODULE m_make_kpoints
