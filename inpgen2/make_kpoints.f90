!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_make_kpoints
  USE m_juDFT
  use m_types_kpts
  USE m_constants
  IMPLICIT NONE
  private
  public :: make_kpoints, add_special_points_default
CONTAINS
  SUBROUTINE make_kpoints(kpts,cell,sym,hybinp,film,l_socorss,bz_integration,l_gamma,str,kptsName,kptsPath)
    USE m_types_kpts
    USE m_types_cell
    USE m_types_sym
    USE m_types_hybinp
    USE m_kpts_kplib
    TYPE(t_kpts),INTENT(out)   :: kpts
    TYPE(t_cell),INTENT(in)    :: cell
    TYPE(t_sym),INTENT(in)     :: sym
    TYPE(t_hybinp), intent(in) :: hybinp
    LOGICAL,INTENT(in)::l_socorss,film
    INTEGER,INTENT(inout)::bz_integration
    LOGICAL,INTENT(inout)::l_gamma
    CHARACTER(len=*),INTENT(inout)::str
    CHARACTER(len=*),INTENT(inout)::kptsName
    CHARACTER(len=*),INTENT(inout)::kptsPath

    LOGICAL:: l_soc_or_ss, l_bzset, l_onlyIdentitySym
    REAL   :: den
    INTEGER:: nk,grid(3)
    character(len=40)::name=""
    !defaults
    l_soc_or_ss=l_socorss

    IF (kptsPath.NE.''.AND.kptsPath.NE.'default') THEN
       CALL set_special_points(kpts,kptsPath)
    ENDIF

    PRINT *,"Processing k-point string: ",str
    !set name
    IF (INDEX(str,"#")>0) THEN
       name=str(:INDEX(str,"#")-1)
       str=trim(adjustl(str(INDEX(str,"#")+1:)))
    ELSE
       name = kptsName
    END IF
    str=ADJUSTL(str)
    DO WHILE(INDEX(str,'@')>0)
       ! Read in the integration method if they are specified on the command line (hist is standard)
       l_bzset = .FALSE. !To warn if there are multiple definitions on the command line
       IF (INDEX(str,'gauss@')==1) THEN
          IF(bz_integration.NE.BZINT_METHOD_HIST) &
              CALL juDFT_warn("You specified the integration method in file and on command line")
          bz_integration=BZINT_METHOD_GAUSS
          IF(l_bzset) &
              CALL juDFT_warn("You specified the integration method multiple times in the command line")
          l_bzset = .TRUE.
          str=str(7:)
       ELSE IF (INDEX(str,'tria@')==1) THEN
          IF(bz_integration.NE.BZINT_METHOD_HIST) &
              CALL juDFT_warn("You specified the integration method in file and on command line")
          bz_integration=BZINT_METHOD_TRIA
          IF(l_bzset) &
              CALL juDFT_warn("You specified the integration method multiple times in the command line")
          l_bzset = .TRUE.
          str=str(6:)
       ELSE IF (INDEX(str,'tetra@')==1) THEN
          IF(bz_integration.NE.BZINT_METHOD_HIST) &
              CALL juDFT_warn("You specified the integration method in file and on command line")
          bz_integration=BZINT_METHOD_TETRA
          IF(l_bzset) &
              CALL juDFT_warn("You specified the integration method multiple times in the command line")
          l_bzset = .TRUE.
          str=str(7:)
       ELSE IF (INDEX(str,'gamma@')==1) THEN
          l_gamma=.TRUE.
          str=str(7:)
       ELSE IF (INDEX(str,'soc@')==1) THEN
          l_soc_or_ss=.TRUE.
          str=str(5:)
       ELSE
          WRITE (*,'(2a)') 'Modifier for k-point set: ', str(1:INDEX(str,'@')-1)
          CALL juDFT_error("Unknown modifier for k-point set.", calledby="make_kpoints")
       ENDIF
    END DO

    l_onlyIdentitySym = .FALSE.
    IF (judft_was_argument("-nosym").OR.judft_was_argument("-noKsym")) THEN
       l_onlyIdentitySym = .TRUE.
    END IF

    l_gamma = l_gamma .or. hybinp%l_hybrid

    IF (INDEX(str,'den=')==1) THEN
       str=str(5:)
       READ(str,*) den
       PRINT *,"Generating a k-point set with density:",den
       CALL init_by_density(kpts,den,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    ELSEIF(INDEX(str,'nk=')==1) THEN
       str=str(4:)
       READ(str,*) nk
       PRINT *,"Generating a k-point set with ",nk," k-points"
       CALL init_by_number(kpts,nk,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    ELSEIF(INDEX(str,'band=')==1) THEN
       str=str(6:)
       READ(str,*) kpts%nkpt
       PRINT *,"Generating a k-point set for bandstructures with ",kpts%nkpt," k-points"
       CALL init_special(kpts,cell,film)
       kpts%kptsKind = KPTS_KIND_PATH
    ELSEIF(INDEX(str,'grid=')==1) THEN
       str=str(6:)
       READ(str,*) grid
       PRINT *,"Generating a k-point grid:",grid
       CALL init_by_grid(kpts,grid,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    ELSEIF(INDEX(str,'file')==1) THEN
       CALL init_by_kptsfile(kpts,film)
       PRINT *,"Reading old kpts file"
    ELSEIF(INDEX(str,'kplib=')==1) THEN
       str=str(7:)
       READ(str,*) kpts%nkpt
       CALL kpts_kplib(cell,sym,kpts)
    ELSEIF(LEN_TRIM(str)<1.OR.INDEX(ADJUSTL(str),'#')==1) THEN
       PRINT *,"Generating default k-point set"
       CALL init_defaults(kpts,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    ELSE
       CALL judft_error(("Could not process -k argument:"//str))
    ENDIF

    if (len_trim(name)>0) kpts%kptsName=name
  END SUBROUTINE make_kpoints

  SUBROUTINE set_special_points(kpts,str)
    IMPLICIT NONE
    TYPE(t_kpts),INTENT(INOUT)::kpts
    CHARACTER(len=*),INTENT(IN)::str

    CHARACTER(len=500) :: rest,l,ll
    INTEGER :: i,err
    real    :: kvec(3)
    rest=str
    DO WHILE(len_TRIM(rest)>1)
       !cut out everything before first ";"
       IF (INDEX(rest,";")>0) THEN
          l=rest(:INDEX(rest,";")-1)
          rest=rest(INDEX(rest,";")+1:)
       ELSE
          l=rest
          rest=""
       ENDIF
       IF (INDEX(l,"=")==0) CALL judft_error("Wrong definition of special k-point:"//l)
       !full definition of special k-point
       ll=l(INDEX(l,"=")+1:)
       READ(ll,*,iostat=err) kvec
       IF (err.NE.0) CALL judft_error("Wrong definition of special k-point:"//l)
       CALL kpts%add_special_line(kvec,l(:INDEX(l,"=")-1))
    END DO
  END SUBROUTINE set_special_points




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

    REAL:: nextp(3),lastp(3)
    REAL,ALLOCATABLE:: d(:)
    INTEGER :: i,ii,iArray(1)
    INTEGER,ALLOCATABLE:: nk(:)
    REAL,ALLOCATABLE :: segmentLengths(:)
    IF (kpts%numSpecialPoints<2) CALL add_special_points_default(kpts,film,cell)
    kpts%nkpt=MAX(kpts%nkpt,kpts%numSpecialPoints)
    !all sepecial kpoints are now set already
    ALLOCATE(nk(kpts%numSpecialPoints-1),d(kpts%numSpecialPoints))
    ALLOCATE(segmentLengths(kpts%numSpecialPoints-1))
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
       nk(i-1)=NINT((kpts%nkpt-kpts%numSpecialPoints)*(d(i)/SUM(d)))
    ENDDO

    DO WHILE (SUM(nk(:))+kpts%numSpecialPoints.NE.kpts%nkpt)
       DO i = 2, kpts%numSpecialPoints
          segmentLengths(i-1) = d(i) / (nk(i-1) + 1)
       END DO
       IF (SUM(nk(:))+kpts%numSpecialPoints.GT.kpts%nkpt) THEN
          iArray = MINLOC(segmentLengths(:))
          nk(iArray(1)) = nk(iArray(1)) - 1
       ELSE
          iArray = MAXLOC(segmentLengths(:))
          nk(iArray(1)) = nk(iArray(1)) + 1
       END IF
    END DO

    ALLOCATE(kpts%bk(3,kpts%numSpecialPoints+SUM(nk)))
    ALLOCATE(kpts%wtkpt(kpts%numSpecialPoints+SUM(nk)))
    kpts%wtkpt = 1.0


    !Generate lines
    kpts%nkpt=1
    DO i=1,kpts%numSpecialPoints-1
       kpts%bk(:,kpts%nkpt)=kpts%specialPoints(:,i)
       kpts%specialPointIndices(i)=kpts%nkpt
       kpts%nkpt=kpts%nkpt+1
       d=(kpts%specialPoints(:,i+1)-kpts%specialPoints(:,i))/(nk(i)+1)
       DO ii=1,nk(i)
          kpts%bk(:,kpts%nkpt)=kpts%specialPoints(:,i)+d*ii
          kpts%nkpt=kpts%nkpt+1
       ENDDO
    ENDDO
    kpts%bk(:,kpts%nkpt)=kpts%specialPoints(:,kpts%numSpecialPoints)
    kpts%specialPointIndices(kpts%numSpecialPoints)=kpts%nkpt
  END SUBROUTINE init_special


  SUBROUTINE init_defaults(kpts,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    USE m_types_cell
    USE m_types_sym
    CLASS(t_kpts),INTENT(out):: kpts
    LOGICAL,INTENT(in)       :: film,l_soc_or_ss,l_gamma
    LOGICAL,INTENT(IN)       :: l_OnlyIdentitySym
    INTEGER,INTENT(in)       :: bz_integration
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym

    INTEGER:: nkpt, nop, nop2

    nop = sym%nop
    nop2 = sym%nop2
    IF(l_OnlyIdentitySym) THEN
       nop = 1
       nop2 = 1
    END IF

    IF (film) THEN
       nkpt = MAX(NINT((3600/cell%area)/nop2),1)
    ELSE
       nkpt = MAX(NINT((216000/cell%omtil)/nop),1)
    ENDIF

    CALL init_by_number(kpts,nkpt,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
  END SUBROUTINE init_defaults

  SUBROUTINE init_by_density(kpts,density,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    USE m_types_cell
    USE m_types_sym
    CLASS(t_kpts),INTENT(out):: kpts
    REAL,INTENT(in)          :: density
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym
    LOGICAL,INTENT(IN)       :: film,l_soc_or_ss,l_gamma
    LOGICAL,INTENT(IN)       :: l_OnlyIdentitySym
    INTEGER,INTENT(IN)       :: bz_integration
    REAL    :: length
    INTEGER :: n,grid(3)

    DO n=1,3
       length=SQRT(DOT_PRODUCT(cell%bmat(n,:),cell%bmat(n,:)))  !TODO why not bmat(:,n)???
       grid(n)=CEILING(density*length)
    END DO
    CALL init_by_grid(kpts,grid,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)

  END SUBROUTINE init_by_density

  SUBROUTINE init_by_number(kpts,nkpt,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    USE m_constants
    USE m_types_cell
    USE m_types_sym
    USE m_types_brZone
    USE m_divi
    USE m_kvecon

    IMPLICIT NONE

    CLASS(t_kpts),INTENT(out):: kpts
    INTEGER,INTENT(IN)       :: nkpt
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym
    LOGICAL,INTENT(IN)       :: film,l_soc_or_ss,l_gamma
    LOGICAL,INTENT(IN)       :: l_OnlyIdentitySym
    INTEGER,INTENT(IN)       :: bz_integration

    TYPE(t_brZone)       :: bz
    INTEGER :: i, j
    INTEGER :: grid(3)

    IF(bz_integration.EQ.BZINT_METHOD_TRIA.AND..NOT.film) THEN
       kpts%nkpt3(:) = 0
       kpts%nkpt = nkpt
       ALLOCATE(kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt))

       WRITE (oUnit,'('' k-points generated with tetrahedron '',''method'')')
       WRITE (oUnit,'(''# k-points generated with tetrahedron '',''method'')')
       WRITE (oUnit,'(3x,'' in irred wedge of 1. Brillouin zone'')')
       WRITE (oUnit,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')

       CALL bz%initBZone(cell, sym, l_soc_or_ss, film,l_OnlyIdentitySym)

       CALL kvecon(kpts%nkpt,mface_const,bz%ncorn,bz%nsym,bz%nface,bz%rltv,bz%fdist,bz%fnorm,bz%cpoint,&
                   kpts%bk)

       DO i = 1, kpts%nkpt
          kpts%wtkpt(i) = 1.0 / kpts%nkpt
          WRITE (oUnit,'(3(f10.7,1x),f12.10,1x,i4,3x,''vkxyz, wghtkp'')') (kpts%bk(j,i),j=1,3),kpts%wtkpt(i),i
       END DO

       DO j=1,kpts%nkpt
          kpts%bk(:,j)=MATMUL(kpts%bk(:,j),cell%amat)/tpi_const
       END DO

       kpts%kptsKind = KPTS_KIND_TRIA_BULK

       RETURN
    END IF
    IF (l_OnlyIdentitySym) THEN
       CALL divi(nkpt,cell%bmat,film,1,1,grid)
    ELSE
       CALL divi(nkpt,cell%bmat,film,sym%nop,sym%nop2,grid)
    END IF
    CALL init_by_grid(kpts,grid,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)

  END SUBROUTINE init_by_number

  SUBROUTINE init_by_grid(kpts,grid,cell,sym,film,bz_integration,l_soc_or_ss,l_gamma,l_OnlyIdentitySym)
    !----------------------------------------------------------------------+
    ! Generate a k-point file with approx. nkpt k-pts or a Monkhorst-Pack  |
    ! set with nmod(i) divisions in i=x,y,z direction. Interface to kptmop |
    ! and kvecon routines of the MD-programm.                              |
    !                                                          G.B. 07/01  |
    !----------------------------------------------------------------------+
    USE m_constants
    USE m_bravais
    USE m_brzone2
    USE m_kptmop
    USE m_kvecon
    USE m_types_cell
    USE m_types_sym
    USE m_types_brZone
    USE m_kptgen_hybrid
    IMPLICIT NONE
    CLASS(t_kpts),INTENT(out):: kpts

    TYPE(t_sym),     INTENT(IN)    :: sym
    TYPE(t_cell),    INTENT(IN)    :: cell
    INTEGER,INTENT(INout)          :: grid(3)
    LOGICAL,INTENT(IN)             :: film,l_soc_or_ss,l_gamma
    LOGICAL,INTENT(IN)             :: l_OnlyIdentitySym
    INTEGER,INTENT(IN)             :: bz_integration

    INTEGER, PARAMETER :: mdir   = 10

    INTEGER ndiv3              ! max. number of tetrahedrons (< 6*(kpts%nkpt+1)

    REAL, ALLOCATABLE    :: vkxyz(:,:)  ! vector of kpoint generated; in cartesian representation
    REAL, ALLOCATABLE    :: wghtkp(:)   !   associated with k-points for BZ integration
    INTEGER, ALLOCATABLE :: ntetra(:,:) ! corners of the tetrahedrons
    REAL, ALLOCATABLE    :: voltet(:)   ! voulmes of the tetrahedrons

    TYPE(t_brZone)       :: bz

    REAL    divis(4)           ! Used to find more accurate representation of k-points

    INTEGER idimens  ! number of dimensions for k-point set (2 or 3)
    INTEGER nreg     ! 1 kpoints in full BZ; 0 kpoints in irrBZ
    INTEGER nfulst   ! 1 kpoints ordered in full stars
    !    (meaningful only for nreg =1; full BZ)
    INTEGER ikzero   ! 0 no shift of k-points;
    ! 1 shift of k-points for better use of sym in irrBZ
    REAL    kzero(3) ! shifting vector to bring one k-point to or
    ! away from (0,0,0) (for even/odd nkpt3)

    INTEGER i,j,k,l,mkpt
    LOGICAL random,l_tria
    REAL as
    REAL binv(3,3)

    kpts%kptsKind = KPTS_KIND_MESH

    IF (l_gamma) THEN
       IF (bz_integration==BZINT_METHOD_TRIA) CALL judft_error("tria and l_gamma incompatible")
       CALL kptgen_hybrid(film,grid,cell,sym,kpts,l_soc_or_ss,l_OnlyIdentitySym)
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


       CALL bz%initBZone(cell, sym, l_soc_or_ss, film,l_OnlyIdentitySym)

       random  = bz_integration==BZINT_METHOD_TRIA.AND..NOT.film

       idimens = MERGE(2,3,film)
       mkpt=PRODUCT(grid(:idimens))

       ALLOCATE (vkxyz(3,mkpt),wghtkp(mkpt))

       IF (bz_integration==BZINT_METHOD_TRIA.AND.random) THEN
          ! Calculate the points for tetrahedron method
          ndiv3 = 6*(mkpt+1)
          ALLOCATE (voltet(ndiv3),ntetra(4,ndiv3))
          kpts%nkpt=mkpt

          WRITE (oUnit,'('' k-points generated with tetrahedron '',''method'')')
          WRITE (oUnit,'(''# k-points generated with tetrahedron '',''method'')')
          WRITE (oUnit,'(3x,'' in irred wedge of 1. Brillouin zone'')')
          WRITE (oUnit,'(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')')

          CALL kvecon(kpts%nkpt,mface_const,bz%ncorn,bz%nsym,bz%nface,bz%rltv,bz%fdist,bz%fnorm,bz%cpoint,&
                      vkxyz)

          DO i = 1, kpts%nkpt
             wghtkp(i) = 1.0 / kpts%nkpt
             WRITE (oUnit,'(3(f10.7,1x),f12.10,1x,i4,3x,''vkxyz, wghtkp'')') (vkxyz(j,i),j=1,3),wghtkp(i),i
          END DO

!          CALL kpttet(kpts%nkpt,ndiv3,&
!               rltv,cell%omtil,nsym,ccr,mdir,mface,&
!               ncorn,nface,fdist,fnorm,cpoint,voltet,ntetra,kpts%ntet,&
!               vkxyz,wghtkp)
       ELSE
          CALL kptmop(bz%idsyst,bz%idtype,grid,&
               bz%rltv,bz%bltv,0,idimens,bz%xvec,bz%fnorm,bz%fdist,bz%ncorn,bz%nface,&
               bz%nedge,bz%cpoint,bz%nsym,bz%ccr,bz%rlsymr,bz%talfa,mkpt,mface_const,mdir,&
               kpts%nkpt,vkxyz,wghtkp)

          WHERE(grid==0) grid=1
       END IF

       DO j=1,kpts%nkpt
          vkxyz(:,j)=MATMUL(vkxyz(:,j),cell%amat)/tpi_const
       END DO

       ALLOCATE(kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt))
       kpts%bk(:,:) = vkxyz(:,:kpts%nkpt)
       kpts%wtkpt(:) = wghtkp(:kpts%nkpt)

    ENDIF

    kpts%nkpt3(:) = grid(:)

  END SUBROUTINE init_by_grid

  SUBROUTINE add_special_points_default(kpts,film,cell,l_check)
    USE m_judft
    USE m_bravais
    USE m_types_cell
    TYPE(t_kpts),INTENT(inout)     :: kpts
    LOGICAL,INTENT(in)             :: film
    LOGICAL,OPTIONAL,INTENT(INOUT) :: l_check
    TYPE(t_cell),INTENT(in)        :: cell

    REAL, PARAMETER :: f12 = 1./2., f14 = 1./4., zro = 0.0
    REAL, PARAMETER :: f34 = 3./4., f38 = 3./8., one = 1.0
    REAL, PARAMETER :: f13 = 1./3., f23 = 2./3.

    INTEGER:: idsyst,idtype
    CALL bravais(cell%amat,idsyst,idtype)

    IF(PRESENT(l_check)) l_check =.FALSE.
    IF (.NOT.film) THEN
       IF ( (idsyst == 1).AND.(idtype ==  3) ) THEN       ! fcc
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/f12,f12,one/) ,"X")
             CALL kpts%add_special_line((/f38,f38,f34/) ,"K")
             CALL kpts%add_special_line((/zro,zro,zro/) ,"g")
             CALL kpts%add_special_line((/f12,f12,f12/) ,"L")
             CALL kpts%add_special_line((/f12,f14,f34/) ,"W")
             CALL kpts%add_special_line((/f12,zro,f12/) ,"X")
             CALL kpts%add_special_line((/zro,zro,zro/) ,"g")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 5).AND.(idtype ==  1) ) THEN       ! rhombohedric (trigonal)
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/f12,f12, f12/) ,"Z")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f14,f14,-f14/) ,"K")
             CALL kpts%add_special_line((/f12,f12,-f12/) ,"Z")
             CALL kpts%add_special_line((/f14,f12,-f14/) ,"W")
             CALL kpts%add_special_line((/zro,f12, zro/) ,"L")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,f12, zro/) ,"F")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 4).AND.(idtype ==  1) ) THEN       ! hexagonal
          IF (cell%bmat(1,1)*cell%bmat(2,1)+cell%bmat(1,2)*cell%bmat(2,2) > 0.0) THEN
             IF(.NOT.PRESENT(l_check)) THEN
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
                CALL kpts%add_special_line((/zro,f12, zro/) ,"M")
                CALL kpts%add_special_line((/f13,f13, zro/) ,"K")
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
                CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
                CALL kpts%add_special_line((/zro,f12, f12/) ,"L")
                CALL kpts%add_special_line((/f13,f13, f12/) ,"H")
                CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
             ELSE
                l_check = .TRUE.
             END IF
          ELSE                                             ! hexagonal (angle = 60)
             IF(.NOT.PRESENT(l_check)) THEN
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
                CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
                CALL kpts%add_special_line((/f13,f23, zro/) ,"K")
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
                CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
                CALL kpts%add_special_line((/f12,f12, f12/) ,"L")
                CALL kpts%add_special_line((/f13,f23, f12/) ,"H")
                CALL kpts%add_special_line((/zro,zro, f12/) ,"A")
             ELSE
                l_check = .TRUE.
             END IF
          ENDIF
       ENDIF
       IF ( (idsyst == 1).AND.(idtype ==  1) ) THEN       ! simple cubic
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,zro, zro/) ,"X")
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
             CALL kpts%add_special_line((/f12,f12, f12/) ,"R")
             CALL kpts%add_special_line((/f12,zro, zro/) ,"X")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,f12, f12/) ,"R")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 1).AND.(idtype ==  2) ) THEN       ! body centered cubic
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,-f12,f12/) ,"H")
             CALL kpts%add_special_line((/zro,zro, f12/) ,"N")
             CALL kpts%add_special_line((/f14,f14, f14/) ,"P")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/zro,zro, f12/) ,"N")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a > c)
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/f12,f12,-f12/) ,"Z")    ! via Lambda and V)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
             CALL kpts%add_special_line((/-f12,f12,f12/) ,"Z")    ! via Y)
             CALL kpts%add_special_line((/zro,zro, f12/) ,"X")    ! via Delta)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/zro,f12, zro/) ,"N")    ! via Q)
             CALL kpts%add_special_line((/f14,f14, f14/) ,"P")    ! via W)
             CALL kpts%add_special_line((/zro,zro, f12/) ,"X")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a < c)
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/-f12,f12,f12/) ,"Z")    ! via F and Sigma)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
             CALL kpts%add_special_line((/zro,zro, f12/) ,"X")    ! via W)
             CALL kpts%add_special_line((/f14,f14, f14/) ,"P")    ! via Q)
             CALL kpts%add_special_line((/zro,f12, zro/) ,"N")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
             CALL kpts%add_special_line((/f12,f12,-f12/) ,"Z")    ! via U and Y)
             CALL kpts%add_special_line((/f12,f12, zro/) ,"X")
             CALL kpts%add_special_line((/f14,f14, f14/) ,"P")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
             CALL kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via Y)
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")    ! via Sigma)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
             CALL kpts%add_special_line((/zro,zro, f12/) ,"Z")    ! via U)
             CALL kpts%add_special_line((/f12,zro, f12/) ,"R")    ! via T)
             CALL kpts%add_special_line((/f12,f12, f12/) ,"A")    ! via S)
             CALL kpts%add_special_line((/zro,zro, f12/) ,"Z")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 3).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          IF(.NOT.PRESENT(l_check)) THEN
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
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
    ELSE
       WRITE(*,*) 'Note:'
       WRITE(*,*) 'Default k point paths for film band structures'
       WRITE(*,*) 'are experimental. If the generated k point path'
       WRITE(*,*) 'is not correct please specify it directly.'
       IF ( (idsyst == 5).AND.(idtype ==  1) ) THEN       ! rhombohedric (trigonal)
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/zro,f12, zro/) ,"L")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,f12, zro/) ,"F")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 4).AND.(idtype ==  1) ) THEN       ! hexagonal
          IF (cell%bmat(1,1)*cell%bmat(2,1)+cell%bmat(1,2)*cell%bmat(2,2) > 0.0) THEN
             IF(.NOT.PRESENT(l_check)) THEN
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
                CALL kpts%add_special_line((/zro,f12, zro/) ,"M")
                CALL kpts%add_special_line((/f13,f13, zro/) ,"K")
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             ELSE
                l_check = .TRUE.
             END IF
          ELSE                                             ! hexagonal (angle = 60)
             IF(.NOT.PRESENT(l_check)) THEN
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
                CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
                CALL kpts%add_special_line((/f13,f23, zro/) ,"K")
                CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             ELSE
                l_check = .TRUE.
             END IF
          ENDIF
       ENDIF
       IF ( (idsyst == 1).AND.(idtype ==  1) ) THEN       ! simple cubic
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")
             CALL kpts%add_special_line((/f12,zro, zro/) ,"X")
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 2).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Delta)
             CALL kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via Y)
             CALL kpts%add_special_line((/f12,f12, zro/) ,"M")    ! via Sigma)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
       IF ( (idsyst == 3).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
          IF(.NOT.PRESENT(l_check)) THEN
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Sigma)
             CALL kpts%add_special_line((/f12,zro, zro/) ,"X")    ! via D)
             CALL kpts%add_special_line((/f12,f12, zro/) ,"S")    ! via C)
             CALL kpts%add_special_line((/zro,f12, zro/) ,"Y")    ! via Delta)
             CALL kpts%add_special_line((/zro,zro, zro/) ,"g")    ! via Lambda)
          ELSE
             l_check = .TRUE.
          END IF
       ENDIF
    END IF

    IF (kpts%numspecialPoints<2) THEN
       IF(.NOT.PRESENT(l_check)) THEN
          CALL judft_error("Not enough special points given and no default found")
       END IF
    END IF
  END SUBROUTINE add_special_points_default
END MODULE m_make_kpoints
