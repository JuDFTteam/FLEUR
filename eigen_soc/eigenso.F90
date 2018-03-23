MODULE m_eigenso
  !
  !*********************************************************************
  !     sets ur and solves the spin-orbit eigenvalue problem in the
  !     second variation procedure.
  !
  !     way: takes e.v. and e.f. from previous scalar-rel. calc.
  !     makes spin-orbit matrix elements solves e.v. and put it on 'eig'
  !
  !     Tree:  eigenso-|- readPotential
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
  SUBROUTINE eigenso(eig_id,mpi,DIMENSION,stars,vacuum,atoms,sphhar,&
                     obsolete,sym,cell,noco,input,kpts,oneD,vTot,enpara)

    USE m_eig66_io, ONLY : read_eig,write_eig
    USE m_spnorb 
    USE m_alineso
    USE m_types
#ifdef CPP_MPI
    USE m_mpi_bc_pot
#endif
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
    TYPE(t_potden),INTENT(IN)    :: vTot
    TYPE(t_enpara),INTENT(IN)    :: enpara
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id       
    !     ..
    !     ..
    !     .. Local Scalars ..
    INTEGER i,j,nk,jspin,n ,l
    INTEGER n_loc,n_plus,i_plus,n_end,nsz,nmat
    LOGICAL l_socvec   !,l_all
    INTEGER wannierspin
    TYPE(t_usdus):: usdus
    !     ..
    !     .. Local Arrays..
    CHARACTER*3 chntype

    REAL,    ALLOCATABLE :: rsopdp(:,:,:,:),rsopdpd(:,:,:,:)
    REAL,    ALLOCATABLE :: rsopp(:,:,:,:),rsoppd(:,:,:,:) 
    REAL,    ALLOCATABLE :: eig_so(:) 
    REAL,    ALLOCATABLE :: rsoplop(:,:,:,:)
    REAL,    ALLOCATABLE :: rsoplopd(:,:,:,:),rsopdplo(:,:,:,:)
    REAL,    ALLOCATABLE :: rsopplo(:,:,:,:),rsoploplop(:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: zso(:,:,:),soangl(:,:,:,:,:,:)

    TYPE(t_mat)::zmat
    TYPE(t_lapw)::lapw

    INTEGER :: ierr
    
    !  ..

    INQUIRE (4649,opened=l_socvec)

    ! To be consistent with angles should be redefined here!
    !noco%theta= -noco%theta
    !noco%phi=   noco%phi+pi_const
    ! now the definition of rotation matrices
    ! is equivalent to the def in the noco-routines

    ALLOCATE(  usdus%us(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd), usdus%dus(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),&
         usdus%uds(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),usdus%duds(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),&
         usdus%ddn(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),&
         usdus%ulos(atoms%nlod,atoms%ntype,DIMENSION%jspd),usdus%dulos(atoms%nlod,atoms%ntype,DIMENSION%jspd),&
         usdus%uulon(atoms%nlod,atoms%ntype,DIMENSION%jspd),usdus%dulon(atoms%nlod,atoms%ntype,DIMENSION%jspd))
   
    IF (input%l_wann.OR.l_socvec) THEN
       wannierspin = 2
    ELSE
       wannierspin = input%jspins
    ENDIF

    !
    !---> set up and solve the eigenvalue problem
    ! --->    radial k-idp s-o matrix elements calc. and storage
    !
#if defined(CPP_MPI)
    !RMA synchronization
    CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
    CALL timestart("eigenso: spnorb")
    !  ..
    ALLOCATE( rsopdp(atoms%ntype,atoms%lmaxd,2,2),rsopdpd(atoms%ntype,atoms%lmaxd,2,2),&
         rsopp(atoms%ntype,atoms%lmaxd,2,2),rsoppd(atoms%ntype,atoms%lmaxd,2,2),&
         rsoplop(atoms%ntype,atoms%nlod,2,2),rsoplopd(atoms%ntype,atoms%nlod,2,2),&
         rsopdplo(atoms%ntype,atoms%nlod,2,2),rsopplo(atoms%ntype,atoms%nlod,2,2),&
         rsoploplop(atoms%ntype,atoms%nlod,atoms%nlod,2,2),&
         soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2) )

    soangl(:,:,:,:,:,:) = CMPLX(0.0,0.0)
    CALL spnorb( atoms,noco,input,mpi, enpara,vTot%mt, rsopp,rsoppd,rsopdp,rsopdpd,usdus,&
         rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop, soangl)
    !
    !Check if SOC is to be scaled for some atom
    DO n=1,atoms%ntype
       IF (ABS(noco%socscale(n)-1.0)>1.E-7) THEN
          IF (mpi%irank==0) WRITE(6,*) "SOC scaled by ",noco%socscale(n)," for atom ",n
          rsopp(n,:,:,:)    =  rsopp(n,:,:,:) * noco%socscale(n)
          rsopdp(n,:,:,:)   =  rsopdp(n,:,:,:)* noco%socscale(n)
          rsoppd(n,:,:,:)   =  rsoppd(n,:,:,:)* noco%socscale(n)
          rsopdpd(n,:,:,:)  =  rsopdpd(n,:,:,:)* noco%socscale(n)
          rsoplop(n,:,:,:)  =  rsoplop(n,:,:,:)* noco%socscale(n)
          rsoplopd(n,:,:,:) =  rsoplopd(n,:,:,:)* noco%socscale(n)
          rsopdplo(n,:,:,:) =  rsopdplo(n,:,:,:)* noco%socscale(n)
          rsopplo(n,:,:,:)  =  rsopplo(n,:,:,:)* noco%socscale(n)
          rsoploplop(n,:,:,:,:) = rsoploplop(n,:,:,:,:)* noco%socscale(n)
       ENDIF
    ENDDO

    IF (mpi%irank==0) THEN
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
    ENDIF
8000 FORMAT (' spin - orbit parameter HR  ')
8001 FORMAT (8f8.4)
9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')
    
 

    IF (mpi%irank==0) THEN
       IF (noco%soc_opt(atoms%ntype+1)) THEN ! .OR. l_all) THEN
!          IF (l_all) THEN
!             WRITE (6,fmt='(A)') 'Only SOC contribution of certain'&
!                  //' atom types included in Hamiltonian.'
!          ELSE 
             WRITE (chntype,'(i3)') atoms%ntype
             WRITE (6,fmt='(A,2x,'//chntype//'l1)') 'SOC contributi'&
                  //'on of certain atom types included in Hamiltonian:',&
                  (noco%soc_opt(n),n=1,atoms%ntype)
!          ENDIF
       ELSE
          WRITE(6,fmt='(A,1x,A)') 'SOC contribution of all atom'//&
               ' types inculded in Hamiltonian.'
       ENDIF
       IF (noco%soc_opt(atoms%ntype+2)) THEN
          WRITE(6,fmt='(A)')&
               'SOC Hamiltonian is constructed by neglecting B_xc.'
       ENDIF
    ENDIF



    ALLOCATE( eig_so(2*DIMENSION%neigd) )
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
       CALL lapw%init(input,noco, kpts,atoms,sym,nk,cell,.FALSE., mpi)
       ALLOCATE( zso(lapw%nv(1)+atoms%nlotot,2*DIMENSION%neigd,wannierspin))
       zso(:,:,:) = CMPLX(0.0,0.0)
       CALL timestart("eigenso: alineso")
       CALL alineso(eig_id,lapw, mpi,DIMENSION,atoms,sym,kpts,&
            input,noco,cell,oneD, rsopp,rsoppd,rsopdp,rsopdpd,nk,&
            rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,&
            usdus,soangl, nsz,nmat, eig_so,zso)
       CALL timestop("eigenso: alineso")
       IF (mpi%irank.EQ.0) THEN
          WRITE (16,FMT=8010) nk,nsz
          WRITE (16,FMT=8020) (eig_so(i),i=1,nsz)
       ENDIF
8010   FORMAT (1x,/,/,' #k=',i6,':',/,&
            ' the',i4,' SOC eigenvalues are:')
8020   FORMAT (5x,5f12.6)

       IF (input%eonly) THEN
          CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz))

       ELSE          
          CALL zmat%alloc(.FALSE.,SIZE(zso,1),nsz)
          DO jspin = 1,wannierspin
             CALL timestart("eigenso: write_eig")  
             zmat%data_c=zso(:,:nsz,jspin)
             CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz),zmat=zmat)

             CALL timestop("eigenso: write_eig")  
          ENDDO
       ENDIF ! (input%eonly) ELSE
       deallocate(zso)
    ENDDO ! DO nk 
    DEALLOCATE (eig_so,rsoploplop,rsopplo,rsopdplo,rsoplopd)
    DEALLOCATE (rsoplop,rsopdp,rsopdpd,rsopp,rsoppd,soangl)


    DEALLOCATE (usdus%us,usdus%dus,usdus%uds,usdus%duds,usdus%ulos,usdus%dulos,usdus%uulon,usdus%dulon,usdus%ddn)
    RETURN
  END SUBROUTINE eigenso
END MODULE m_eigenso
