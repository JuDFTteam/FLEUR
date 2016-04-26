MODULE m_eigenso
  !
  !*********************************************************************
  !     sets up and solves the spin-orbit eigenvalue problem in the
  !     second variation procedure.
  !
  !     way: takes e.v. and e.f. from previous scalar-rel. calc.
  !     makes spin-orbit matrix elements solves e.v. and put it on 'eig'
  !
  !     Tree:  eigenso-|- loddop
  !                    |- spnorb  : sets up s-o parameters 
  !                    |    |- soinit - sorad  : radial part
  !                    |    |- sgml            : diagonal angular parts
  !                    |    |- anglso          : non-diagonal -"-
  !                    |
  !                    |- alineso : sets up and solves e.v. problem
  !                         |- hsohelp
  !                         |- hsoham
  !
  !**********************************************************************
  !
CONTAINS
  SUBROUTINE eigenso(eig_id, mpi,DIMENSION,stars,vacuum,atoms,sphhar,&
       obsolete,sym,cell,noco, input,kpts, oneD)
    !
    USE m_eig66_io, ONLY : read_eig,write_eig
    USE m_spnorb 
    USE m_alineso
    USE m_loddop
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_obsolete),INTENT(IN)  :: obsolete
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id       
    !     ..
    !     ..
    !     .. Local Scalars ..
    INTEGER i,j,nk,jspin ,iter ,n ,l
    INTEGER n_loc,n_plus,i_plus,n_end,nsz
    LOGICAL l_all,l_file,l_socvec
    INTEGER wannierspin
    TYPE(t_enpara) :: enpara
    TYPE(t_usdus):: usdus
    !     ..
    !     .. Local Arrays..
    CHARACTER*3 chntype

    INTEGER, ALLOCATABLE :: kveclo(:)
    REAL,    ALLOCATABLE :: rsopdp(:,:,:,:),rsopdpd(:,:,:,:)
    REAL,    ALLOCATABLE :: rsopp(:,:,:,:),rsoppd(:,:,:,:) 
    REAL,    ALLOCATABLE :: eig_so(:) 
    REAL,    ALLOCATABLE :: rsoplop(:,:,:,:)
    REAL,    ALLOCATABLE :: rsoplopd(:,:,:,:),rsopdplo(:,:,:,:)
    REAL,    ALLOCATABLE :: rsopplo(:,:,:,:),rsoploplop(:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: zso(:,:,:),soangl(:,:,:,:,:,:)

    REAL,    ALLOCATABLE :: vz(:,:,:),vr(:,:,:,:)
    COMPLEX, ALLOCATABLE :: vzxy(:,:,:,:),vpw(:,:)



    !  ..

    INQUIRE (4649,opened=l_socvec)

    ! To be consistent with angles should be redefined here!
    !noco%theta= -noco%theta
    !noco%phi=   noco%phi+pi_const
    ! now the definition of rotation matrices
    ! is equivalent to the def in the noco-routines
    !
    ! load potential from file pottot (=unit 8)
    !
    ALLOCATE ( vz(vacuum%nmzd,2,DIMENSION%jspd),vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,DIMENSION%jspd),&
         vzxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,DIMENSION%jspd),vpw(stars%n3d,DIMENSION%jspd) )

    OPEN (8,file='pottot',form='unformatted',status='old')
    CALL loddop(&
         stars,vacuum,atoms,sphhar,&
         input,sym,&
         8,&
         iter,vr,vpw,vz,vzxy)
    CLOSE(8)

    DEALLOCATE ( vz,vzxy,vpw )

    ALLOCATE(  usdus%us(0:atoms%lmaxd,atoms%ntypd,DIMENSION%jspd), usdus%dus(0:atoms%lmaxd,atoms%ntypd,DIMENSION%jspd),&
         usdus%uds(0:atoms%lmaxd,atoms%ntypd,DIMENSION%jspd),usdus%duds(0:atoms%lmaxd,atoms%ntypd,DIMENSION%jspd),&
         usdus%ddn(0:atoms%lmaxd,atoms%ntypd,DIMENSION%jspd),kveclo(atoms%nlotot),&
         usdus%ulos(atoms%nlod,atoms%ntypd,DIMENSION%jspd),usdus%dulos(atoms%nlod,atoms%ntypd,DIMENSION%jspd),&
         usdus%uulon(atoms%nlod,atoms%ntypd,DIMENSION%jspd),usdus%dulon(atoms%nlod,atoms%ntypd,DIMENSION%jspd),&
         enpara%evac0(2,DIMENSION%jspd),enpara%ello0(atoms%nlod,atoms%ntypd,DIMENSION%jspd),&
         enpara%el0(0:atoms%lmaxd,atoms%ntypd,DIMENSION%jspd))

    INQUIRE (file='wann_inp',exist=l_file)
    IF (l_file.OR.l_socvec) THEN
       wannierspin = 2
    ELSE
       wannierspin = input%jspins
    ENDIF

    !
    !---> set up and solve the eigenvalue problem
    !
    !--->    radial k-idp s-o matrix elements calc. and storage
    !
    DO jspin = 1, input%jspins
       CALL read_eig(eig_id,&
            1,jspin,&
            el=enpara%el0(:,:,jspin),&
            ello=enpara%ello0(:,:,jspin),evac=enpara%evac0(:,jspin))
    ENDDO
    CALL timestart("eigenso: spnorb")
    !  ..
    ALLOCATE( rsopdp(atoms%ntypd,atoms%lmaxd,2,2),rsopdpd(atoms%ntypd,atoms%lmaxd,2,2),&
         rsopp(atoms%ntypd,atoms%lmaxd,2,2),rsoppd(atoms%ntypd,atoms%lmaxd,2,2),&
         rsoplop(atoms%ntypd,atoms%nlod,2,2),rsoplopd(atoms%ntypd,atoms%nlod,2,2),&
         rsopdplo(atoms%ntypd,atoms%nlod,2,2),rsopplo(atoms%ntypd,atoms%nlod,2,2),&
         rsoploplop(atoms%ntypd,atoms%nlod,atoms%nlod,2,2),&
         soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2) )

    soangl(:,:,:,:,:,:) = CMPLX(0.0,0.0)
    CALL spnorb(&
         atoms,noco,input,mpi,&
         enpara,vr,&
         input%sso_opt(atoms%ntype+2), &
         rsopp,rsoppd,rsopdp,rsopdpd,usdus,&
         rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,&
         soangl)
    !
    l_all = .FALSE.
    INQUIRE (file='allbut',exist=l_all)
    IF (l_all) THEN
       OPEN (1,file='allbut',form='formatted')
       READ (1,*) n
       WRITE (*,*) 'allbut',n
       CLOSE (1)
       rsopp(1:n-1,:,:,:) = 0.0 ; rsopp(n+1:atoms%ntypd,:,:,:) = 0.0 
       rsopdp(1:n-1,:,:,:) = 0.0 ; rsopdp(n+1:atoms%ntypd,:,:,:) = 0.0 
       rsoppd(1:n-1,:,:,:) = 0.0 ; rsoppd(n+1:atoms%ntypd,:,:,:) = 0.0 
       rsopdpd(1:n-1,:,:,:) = 0.0 ; rsopdpd(n+1:atoms%ntypd,:,:,:) = 0.0 
       rsoplop(1:n-1,:,:,:) = 0.0 ; rsoplop(n+1:atoms%ntypd,:,:,:) = 0.0 
       rsoplopd(1:n-1,:,:,:) = 0.0 ; rsoplopd(n+1:atoms%ntypd,:,:,:) = 0.0 
       rsopdplo(1:n-1,:,:,:) = 0.0 ; rsopdplo(n+1:atoms%ntypd,:,:,:) = 0.0 
       rsopplo(1:n-1,:,:,:) = 0.0 ; rsopplo(n+1:atoms%ntypd,:,:,:) = 0.0
       rsoploplop(1:n-1,:,:,:,:) = 0.0 ; rsoploplop(n+1:atoms%ntypd,:,:,:,:) = 0.0
    ENDIF
    l_all = .FALSE.
    INQUIRE (file='socscale',exist=l_all)
    IF (l_all) THEN
       OPEN (1,file='socsacle',form='formatted')
       READ (1,*) n
       WRITE (*,*) 'SOC scaled by ',n,"%"
       CLOSE (1)
       rsopp(:,:,:,:) = n/100.* rsopp
       rsopdp(:,:,:,:) =  n/100.*rsopdp
       rsoppd(:,:,:,:) =  n/100.*rsoppd
       rsopdpd(:,:,:,:) =  n/100.*rsopdpd
       rsoplop(:,:,:,:) =  n/100.*rsoplop
       rsoplopd(:,:,:,:) =  n/100.*rsoplopd
       rsopdplo(:,:,:,:) =  n/100.*rsopdplo
       rsopplo(:,:,:,:) =  n/100.* rsopplo
       rsoploplop(:,:,:,:,:) = n/100.*rsoploplop

       DO n = 1,atoms%ntype
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopp(n,l,2,1),l=1,3)
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsoppd(n,l,2,1),l=1,3)
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdp(n,l,2,1),l=1,3)
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsopdpd(n,l,2,1),l=1,3)
       ENDDO
8000   FORMAT (' spin - orbit parameter HR  ')
8001   FORMAT (8f8.4)
9000   FORMAT (5x,' p ',5x,' d ', 5x, ' f ')

    ENDIF


    IF (mpi%irank==0) THEN
       IF (noco%soc_opt(atoms%ntype+1) .OR. l_all) THEN
          IF (l_all) THEN
             WRITE (6,fmt='(A)') 'Only SOC contribution of certain'&
                  //' atom types included in Hamiltonian.'
          ELSE 
             WRITE (chntype,'(i3)') atoms%ntype
             WRITE (6,fmt='(A,2x,'//chntype//'l1)') 'SOC contributi'&
                  //'on of certain atom types included in Hamiltonian:',&
                  (noco%soc_opt(n),n=1,atoms%ntype)
          ENDIF
       ELSE
          WRITE(6,fmt='(A,1x,A)') 'SOC contribution of all atom'//&
               ' types inculded in Hamiltonian.'
       ENDIF
       IF (noco%soc_opt(atoms%ntype+2)) THEN
          WRITE(6,fmt='(A)')&
               'SOC Hamiltonian is constructed by neglecting B_xc.'
       ENDIF
    ENDIF



    ALLOCATE( zso(DIMENSION%nbasfcn,2*DIMENSION%neigd,wannierspin),eig_so(2*DIMENSION%neigd) )
    zso(:,:,:) = CMPLX(0.0,0.0)
    soangl(:,:,:,:,:,:) = CONJG(soangl(:,:,:,:,:,:))
    CALL timestop("eigenso: spnorb")
    !
    !--->    loop over k-points: each can be a separate task
    !
    n_loc = INT(kpts%nkpt/mpi%isize)
    n_plus = kpts%nkpt - mpi%isize*n_loc
    i_plus = -1
    IF (mpi%irank.LT.n_plus) i_plus = 0
    n_end = (mpi%irank+1)+(n_loc+i_plus)*mpi%isize
    !
    !--->  start loop k-pts
    !
    DO  nk = mpi%irank+1,n_end,mpi%isize

       CALL timestart("eigenso: alineso")
       CALL alineso(eig_id,&
            mpi,DIMENSION,atoms,sym,&
            input,noco,cell,oneD,&
            rsopp,rsoppd,rsopdp,rsopdpd,nk,&
            rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,&
            usdus,soangl,&
            kveclo,enpara%ello0,nsz,&
            eig_so,zso)
       CALL timestop("eigenso: alineso")
       IF (mpi%irank.EQ.0) THEN
          WRITE (16,FMT=8010) nk,nsz
          WRITE (16,FMT=8020) (eig_so(i),i=1,nsz)
       ENDIF
8010   FORMAT (1x,/,/,' #k=',i6,':',/,&
            ' the',i4,' SOC eigenvalues are:')
8020   FORMAT (5x,5f12.6)

       IF (input%eonly) THEN
          CALL write_eig(eig_id,&
               nk,jspin,neig=nsz,neig_total=nsz,nmat=SIZE(zso,1),&
               eig=eig_so(:nsz))

       ELSE

          DO jspin = 1,wannierspin
             CALL write_eig(eig_id,&
                  nk,jspin,neig=nsz,neig_total=nsz,nmat=SIZE(zso,1),&
                  eig=eig_so(:nsz),z=zso(:,:nsz,jspin))

             CALL timestop("eigenso: write_eig")  
          ENDDO

       ENDIF ! (input%eonly) ELSE

    ENDDO ! DO nk 
    DEALLOCATE (zso,eig_so,rsoploplop,rsopplo,rsopdplo,rsoplopd)
    DEALLOCATE (rsoplop,rsopdp,rsopdpd,rsopp,rsoppd,soangl)


    DEALLOCATE ( vr,usdus%us,usdus%dus,usdus%uds,usdus%duds,usdus%ulos,usdus%dulos,usdus%uulon,usdus%dulon,usdus%ddn )
    RETURN
  END SUBROUTINE eigenso
END MODULE m_eigenso
