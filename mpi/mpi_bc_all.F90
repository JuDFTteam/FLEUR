MODULE m_mpi_bc_all
CONTAINS
  SUBROUTINE mpi_bc_all(&
       mpi,stars,sphhar,atoms,obsolete,sym,&
       kpts,jij,dimension,input,banddos,sliceplot,&
       vacuum,cell,enpara,noco,oneD,&
       MM,nn2d_1,n2d_1,nop_1,&
        xcpot,hybrid)
    !
    !**********************************************************************
    USE m_hybridmix
    USE m_types
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    TYPE(t_xcpot),INTENT(INOUT)      :: xcpot
    TYPE(t_mpi),INTENT(INOUT)        :: mpi
    TYPE(t_dimension),INTENT(INOUT)  :: dimension
    TYPE(t_oneD),INTENT(INOUT)       :: oneD
    TYPE(t_hybrid),INTENT(INOUT)     :: hybrid
    TYPE(t_enpara),INTENT(INOUT)     :: enpara
    TYPE(t_obsolete),INTENT(INOUT)   :: obsolete
    TYPE(t_banddos),INTENT(INOUT)    :: banddos
    TYPE(t_sliceplot),INTENT(INOUT)  :: sliceplot
    TYPE(t_input),INTENT(INOUT)      :: input
    TYPE(t_vacuum),INTENT(INOUT)     :: vacuum
    TYPE(t_noco),INTENT(INOUT)       :: noco
    TYPE(t_jij),INTENT(INOUT)        :: jij
    TYPE(t_sym),INTENT(INOUT)        :: sym
    TYPE(t_stars),INTENT(INOUT)      :: stars
    TYPE(t_cell),INTENT(INOUT)       :: cell
    TYPE(t_kpts),INTENT(INOUT)       :: kpts
    TYPE(t_sphhar),INTENT(INOUT)     :: sphhar
    TYPE(t_atoms),INTENT(INOUT)      :: atoms
    !     .. Scalar Arguments ..
    INTEGER,INTENT(INOUT) :: n2d_1,nn2d_1,nop_1,MM
    !     ..
    !     .. Local Scalars ..
    INTEGER n
    REAL rdum
    !     .. Local Arrays ..
    INTEGER i(36),ierr(3)
    REAL    r(28)
    LOGICAL l(43)
    !     ..
    !     .. External Subroutines.. 
    EXTERNAL MPI_BCAST


    mpi%typ_real=MPI_DOUBLE_PRECISION
    mpi%typ_complex=MPI_DOUBLE_COMPLEX

    IF (mpi%irank.EQ.0) THEN
       i(1)=1 ; i(2)=obsolete%lpr ; i(3)=atoms%ntype ; i(5)=obsolete%nwd ; i(6)=input%isec1
       i(7)=stars%ng2 ; i(8)=stars%nq3 ; i(9)=vacuum%nmz ; i(10)=vacuum%nmzxy ; i(11)=obsolete%lepr 
       i(12)=input%jspins ; i(13)=vacuum%nvac ; i(14)=input%itmax ; i(15)=sliceplot%kk ; i(16)=vacuum%layers
       i(17)=sliceplot%nnne ; i(18)=banddos%ndir ; i(19)=stars%mx1 ; i(20)=stars%mx2 ; i(21)=stars%mx3
       i(22)=atoms%n_u ; i(23) = sym%nop2 ; i(24) = sym%nsymt ; i(25) = xcpot%icorr
       i(26)=vacuum%nstars ; i(27)=vacuum%nstm ; i(28)=oneD%odd%nq2 ; i(29)=oneD%odd%nop
       i(30)=input%gw ; i(31)=input%gw_neigd ; i(32)=hybrid%ewaldlambda ; i(33)=hybrid%lexp 
       i(34)=hybrid%bands1 ; i(35)=hybrid%bands2 ; i(36)=input%imix
       r(1)=cell%omtil ; r(2)=cell%area ; r(3)=vacuum%delz ; r(4)=cell%z1 ; r(5)=input%alpha
       r(6)=sliceplot%e1s ; r(7)=sliceplot%e2s ; r(8)=noco%theta ; r(9)=noco%phi ; r(10)=vacuum%tworkf 
       r(11)=vacuum%locx(1) ; r(12)=vacuum%locx(2); r(13)=vacuum%locy(1) ; r(14)=vacuum%locy(2)
       r(15)=input%efield%sigma ; r(16)=input%efield%zsigma ; r(17)=noco%mix_b; r(18)=cell%vol
       r(19)=cell%volint ; r(20)=hybrid%gcutm1 ; r(21)=hybrid%tolerance1 ; r(22)=hybrid%gcutm2
       r(23)=hybrid%tolerance2 ; r(24)=input%delgau ; r(25)=input%tkb ; r(26)=input%efield%vslope
       r(27)=aMix_VHSE() ; r(28)=omega_VHSE()
       l(1)=input%eonly ; l(2)=obsolete%form66 ; l(3)=input%secvar ; l(4)=sym%zrfs ; l(5)=input%film
       l(6)=sym%invs ; l(7)=sym%invs2 ; l(8)=input%l_bmt ; l(9)=input%l_f ; l(10)=input%cdinf
       l(11)=banddos%dos ;l(12)=sym%l_zref ; l(13)=banddos%vacdos ; l(14)=input%integ ; l(15)=sliceplot%iplot
       l(16)=input%strho ; l(17)=input%swsp ; l(18)=input%lflip ; l(19)=obsolete%l_f2u ; l(20)=obsolete%l_u2f
       l(21)=input%pallst ; l(22)=sliceplot%slice ; l(23)=noco%l_soc ; l(24)=vacuum%starcoeff
       l(25)=noco%l_noco ; l(26)=noco%l_ss; l(27)=noco%l_mperp; l(28)=noco%l_constr
       l(29)=oneD%odd%d1 ; l(30)=jij%l_J ; l(31)=jij%l_disp ; l(32)=input%ctail
       l(33)=obsolete%eig66(1) ; l(34)=obsolete%eig66(2) ; l(35)=input%sso_opt(1)
       l(36)=input%sso_opt(2) ; l(37)=obsolete%pot8; l(38)=input%efield%l_segmented
       l(39)=sym%symor ; l(40)=input%frcor ; l(41)=input%tria ; l(42)=input%efield%dirichlet
       l(43)=input%efield%l_dirichlet_coeff
    ENDIF
    !
    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,0,mpi%mpi_comm,ierr)
    hybrid%bands1=i(34) ; hybrid%bands2=i(35) ; input%imix=i(36)
    input%gw=i(30) ; input%gw_neigd=i(31) ; hybrid%ewaldlambda=i(32) ; hybrid%lexp=i(33)
    vacuum%nstars=i(26) ; vacuum%nstm=i(27) ; oneD%odd%nq2=i(28) ; oneD%odd%nop=i(29)
    atoms%n_u=i(22) ; sym%nop2=i(23) ; sym%nsymt = i(24) ; xcpot%icorr=i(25)
    sliceplot%nnne=i(17) ; banddos%ndir=i(18) ; stars%mx1=i(19) ; stars%mx2=i(20) ; stars%mx3=i(21)
    input%jspins=i(12) ; vacuum%nvac=i(13) ; input%itmax=i(14) ; sliceplot%kk=i(15) ; vacuum%layers=i(16)
    stars%ng2=i(7) ; stars%nq3=i(8) ; vacuum%nmz=i(9) ; vacuum%nmzxy=i(10) ; obsolete%lepr=i(11)
    obsolete%lpr=i(2) ; atoms%ntype=i(3) ; obsolete%nwd=i(5) ; input%isec1=i(6)
    !
    CALL MPI_BCAST(r,SIZE(r),mpi%type_real,0,mpi%mpi_comm,ierr)
    rdum=aMix_VHSE( r(27) ); rdum=omega_VHSE( r(28) )
    hybrid%tolerance2=r(23) ; input%delgau=r(24) ; input%tkb=r(25) ; input%efield%vslope=r(26)
    cell%volint=r(19) ; hybrid%gcutm1=r(20) ; hybrid%tolerance1=r(21) ; hybrid%gcutm2=r(22)
    input%efield%sigma=r(15) ; input%efield%zsigma=r(16); noco%mix_b=r(17); cell%vol=r(18);
    vacuum%locx(1)=r(11); vacuum%locx(2)=r(12); vacuum%locy(1)=r(13); vacuum%locy(2)=r(14)
    sliceplot%e1s=r(6) ; sliceplot%e2s=r(7) ; noco%theta=r(8) ; noco%phi=r(9) ; vacuum%tworkf=r(10)
    cell%omtil=r(1) ; cell%area=r(2) ; vacuum%delz=r(3) ; cell%z1=r(4) ; input%alpha=r(5)
    !
    CALL MPI_BCAST(l,SIZE(l),MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    input%efield%l_dirichlet_coeff = l(43)
    sym%symor=l(39) ; input%frcor=l(40) ; input%tria=l(41) ; input%efield%dirichlet = l(42)
    input%sso_opt(2)=l(36) ; obsolete%pot8=l(37) ; input%efield%l_segmented=l(38)
    obsolete%eig66(1)=l(33) ; obsolete%eig66(2)=l(34) ; input%sso_opt(1)=l(35)
    oneD%odd%d1=l(29) ; jij%l_J=l(30) ; jij%l_disp=l(31) ; input%ctail=l(32)
    noco%l_noco=l(25) ; noco%l_ss=l(26) ; noco%l_mperp=l(27) ; noco%l_constr=l(28)
    input%pallst=l(21) ; sliceplot%slice=l(22) ; noco%l_soc=l(23) ; vacuum%starcoeff=l(24)
    input%strho=l(16) ; input%swsp=l(17) ; input%lflip=l(18) ; obsolete%l_f2u=l(19) ; obsolete%l_u2f=l(20)
    banddos%dos=l(11) ;sym%l_zref=l(12) ; banddos%vacdos=l(13) ; input%integ=l(14) ; sliceplot%iplot=l(15)
    sym%invs=l(6) ; sym%invs2=l(7) ; input%l_bmt=l(8) ; input%l_f=l(9) ; input%cdinf=l(10)
    input%eonly=l(1) ; obsolete%form66=l(2) ; input%secvar=l(3) ; sym%zrfs=l(4) ; input%film=l(5)
    input%efield%l_segmented = l(38) ; sym%symor=l(39); input%efield%dirichlet = l(40)
    input%efield%l_dirichlet_coeff = l(41)
    !
    ! -> Broadcast the arrays:
    IF (input%efield%l_segmented) THEN
       IF (.NOT. ALLOCATED (input%efield%rhoEF))&
            &    ALLOCATE (input%efield%rhoEF(3*stars%k1d*3*stars%k2d-1,vacuum%nvac))
       n = (3*stars%k1d*3*stars%k2d-1)*vacuum%nvac
       CALL MPI_BCAST (input%efield%rhoEF,n,MPI_REAL,0,mpi%mpi_comm,ierr)
    END IF
    IF (input%efield%l_dirichlet_coeff) THEN
       IF (.NOT. ALLOCATED (input%efield%C1)) THEN
          ALLOCATE (input%efield%C1(stars%ng2-1))
          ALLOCATE (input%efield%C2(stars%ng2-1))
       END IF
       n = stars%ng2-1
       CALL MPI_BCAST (input%efield%C1,n,MPI_REAL,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST (input%efield%C2,n,MPI_REAL,0,mpi%mpi_comm,ierr)
    END IF

    CALL MPI_BCAST(stars%ustep,stars%n3d,mpi%typ_complex,0,mpi%mpi_comm,ierr)
    n = sphhar%memd*(sphhar%nlhd+1)*sphhar%ntypsd
    CALL MPI_BCAST(sphhar%clnu,n,mpi%typ_complex,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sphhar%mlh,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sphhar%nlh,sphhar%ntypsd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n = (sphhar%nlhd+1)*sphhar%ntypsd
    CALL MPI_BCAST(sphhar%nmem,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sphhar%llh,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%jri,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ncv,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ntypsy,atoms%natd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%neq,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lnonsph,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lmax,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%invsat,atoms%natd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invsatnr,atoms%natd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ngopr,atoms%natd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%mrot,9*sym%nop,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%ig2,stars%n3d,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n = (2*stars%k1d+1)*(2*stars%k2d+1)*(2*stars%k3d+1)
    CALL MPI_BCAST(stars%ig,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%rgphs,n,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%ellow,obsolete%nwdd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%elup,obsolete%nwdd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%rkmax,obsolete%nwdd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%rmt,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%volmts,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%dx,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%sk3,stars%n3d,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(enpara%evac0,2*dimension%jspd*obsolete%nwdd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(cell%amat,9,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(cell%bmat,9,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(cell%bbmat,9,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%taual,3*atoms%natd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%pos,3*atoms%natd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%tau,3*sym%nop,mpi%type_real,0,mpi%mpi_comm,ierr)
    n = (atoms%lmaxd+1)*atoms%ntypd*dimension%jspd*obsolete%nwdd
    CALL MPI_BCAST(enpara%el0,n,mpi%type_real,0,mpi%mpi_comm,ierr)
    n = atoms%nlod*atoms%ntypd*dimension%jspd
    CALL MPI_BCAST(enpara%ello0,n,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%rmsh,atoms%jmtd*atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    !
    CALL MPI_BCAST(kpts%nkpt,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%bk,3*kpts%nkptd*obsolete%nwdd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%wtkpt,kpts%nkptd*obsolete%nwdd,mpi%type_real,0,mpi%mpi_comm,ierr)
    !
    n = atoms%natd*sym%nop
    CALL MPI_BCAST(sym%invarop,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%multab,sym%nop**2,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invarind,atoms%natd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invtab,sym%nop,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invsatnr,atoms%natd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%kq2_fft,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%kq3_fft,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%kv2,2*stars%n2d,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%kv3,3*stars%n3d,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%ng3_fft,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%ng2_fft,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%kmxq_fft,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%kmxq2_fft,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(vacuum%izlay,vacuum%layerd*2,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%nstr,stars%n3d,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%nstr2,stars%n2d,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%igfft,dimension%nn3d*2,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%kq1_fft,obsolete%nwdd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%pgfft,dimension%nn3d,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%zatom,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%efield%sig_b,2,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%zelec,obsolete%nwdd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ncst,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%nlo,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n =  atoms%nlod*atoms%ntypd
    CALL MPI_BCAST(atoms%llo,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ulo_der,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%l_dulo,n,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    n = (atoms%llod+1)*atoms%ntypd
    CALL MPI_BCAST(atoms%lo1l,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%nlol,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n = dimension%jspd*atoms%ntypd
    CALL MPI_BCAST(enpara%skiplo,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%alph,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%beta,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%b_con,atoms%ntypd*2,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%l_relax,atoms%ntypd,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%l_geo,atoms%ntypd,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%soc_opt,atoms%ntypd+2,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%qss,3,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%l,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%u,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%j,atoms%ntypd,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms(:)%l_amf,MPI_LOGICAL,0,&
         &                                           mpi,ierr)
    CALL MPI_BCAST(atoms%lapw_l,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n = 7*7*3*sym%nop
    CALL MPI_BCAST(sym%d_wgn,n,mpi%typ_complex,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(oneD%nstr1,n2d_1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(stars%sk2,stars%n2d,mpi%type_real,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(oneD%tau1,3*nop_1,mpi%type_real,0,mpi%mpi_comm,ierr)
    IF (oneD%odd%d1) THEN
       CALL MPI_BCAST(stars%phi2,stars%n2d,mpi%type_real,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%tau1,3*nop_1,mpi%type_real,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%mrot1,9*nop_1,mpi%type_real,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%kv1,2*n2d_1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%ngopr1,atoms%natd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%igfft1,nn2d_1*2,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%pgfft1,nn2d_1,mpi%type_real,0,mpi%mpi_comm,ierr)
       n = (2*stars%k3d + 1)*(2*MM +1)
       CALL MPI_BCAST(oneD%ig1,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%invtab1,nop_1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%multab1,2*nop_1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    ENDIF
    !--- J<
    IF (jij%l_J) THEN
       CALL MPI_BCAST(jij%nmagn,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(jij%mtypes,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(jij%thetaJ,1,mpi%type_real,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(jij%qj,3*jij%nqptd,mpi%type_real,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(jij%l_magn,atoms%ntypd,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(jij%magtype,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(jij%nmagtype,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    ENDIF
    !--- J>
    !--- HF<
    CALL MPI_BCAST(hybrid%lcutwf,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%nkpt3,3,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(hybrid%select1,4*atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(hybrid%lcutm1,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(hybrid%select2,4*atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(hybrid%lcutm2,atoms%ntypd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    !--- HF>
    !
    !

    RETURN
  END SUBROUTINE mpi_bc_all
END MODULE m_mpi_bc_all
