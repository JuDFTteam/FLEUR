!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_all
CONTAINS
  SUBROUTINE mpi_bc_all(&
       mpi,stars,sphhar,atoms,obsolete,sym,&
       kpts,DIMENSION,input,field,banddos,sliceplot,&
       vacuum,cell,enpara,noco,oneD,&
        hybrid)
    !
    !**********************************************************************
    USE m_types
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    TYPE(t_mpi),INTENT(INOUT)        :: mpi
    TYPE(t_dimension),INTENT(INOUT)  :: dimension
    TYPE(t_oneD),INTENT(INOUT)       :: oneD
    TYPE(t_hybrid),INTENT(INOUT)     :: hybrid
    TYPE(t_enpara),INTENT(INOUT)     :: enpara
    TYPE(t_obsolete),INTENT(INOUT)   :: obsolete
    TYPE(t_banddos),INTENT(INOUT)    :: banddos
    TYPE(t_sliceplot),INTENT(INOUT)  :: sliceplot
    TYPE(t_input),INTENT(INOUT)      :: input
    TYPE(t_field),INTENT(INOUT)      :: field
    TYPE(t_vacuum),INTENT(INOUT)     :: vacuum
    TYPE(t_noco),INTENT(INOUT)       :: noco
    TYPE(t_sym),INTENT(INOUT)        :: sym
    TYPE(t_stars),INTENT(INOUT)      :: stars
    TYPE(t_cell),INTENT(INOUT)       :: cell
    TYPE(t_kpts),INTENT(INOUT)       :: kpts
    TYPE(t_sphhar),INTENT(INOUT)     :: sphhar
    TYPE(t_atoms),INTENT(INOUT)      :: atoms
    !     .. Scalar Arguments ..
    INTEGER n
    REAL rdum
    !     .. Local Arrays ..
    INTEGER i(43),ierr(3)
    REAL    r(34)
    LOGICAL l(46)
    !     ..
    !     .. External Subroutines..
#ifdef CPP_MPI    
    EXTERNAL MPI_BCAST

    call priv_mpi_bc_stars(mpi,stars)
    
    IF (mpi%irank.EQ.0) THEN
       i(1)=1 ; i(2)=input%coretail_lmax;i(3)=atoms%ntype  ; i(5)=1
       i(7)=stars%ng2 ; i(8)=stars%ng3 ; i(9)=vacuum%nmz ; i(10)=vacuum%nmzxy ; i(11)=obsolete%lepr 
       i(12)=input%jspins ; i(13)=vacuum%nvac ; i(14)=input%itmax ; i(15)=sliceplot%kk ; i(16)=vacuum%layers
       i(17)=sliceplot%nnne ; i(18)=banddos%ndir ; i(19)=stars%mx1 ; i(20)=stars%mx2 ; i(21)=stars%mx3
       i(22)=atoms%n_u ; i(23) = sym%nop2 ; i(24) = sym%nsymt ; i(25) = stars%kimax ; i(26) = stars%kimax2
       i(27)=vacuum%nstars ; i(28)=vacuum%nstm ; i(29)=oneD%odd%nq2 ; i(30)=oneD%odd%nop
       i(31)=input%gw ; i(32)=input%gw_neigd ; i(33)=hybrid%ewaldlambda ; i(34)=hybrid%lexp 
       i(35)=hybrid%bands1 ; i(36)=input%maxiter ; i(37)=input%imix ; i(38)=banddos%orbCompAtom
       i(39)=input%kcrel;i(40)=banddos%s_cell_x;i(41)=banddos%s_cell_y;i(42)=banddos%s_cell_z; i(43)=sliceplot%iplot

       r(1)=cell%omtil ; r(2)=cell%area ; r(3)=vacuum%delz ; r(4)=cell%z1 ; r(5)=input%alpha
       r(6)=sliceplot%e1s ; r(7)=sliceplot%e2s ; r(8)=noco%theta; r(9)=noco%phi; r(10)=vacuum%tworkf 
       r(11)=vacuum%locx(1) ; r(12)=vacuum%locx(2); r(13)=vacuum%locy(1) ; r(14)=vacuum%locy(2)
       r(15)=input%sigma ; r(16)=field%efield%zsigma ; r(17)=noco%mix_b; r(18)=cell%vol
       r(19)=cell%volint ; r(20)=hybrid%gcutm1 ; r(21)=hybrid%tolerance1 ; r(22)=0.0
       r(23)=0.0 ; r(24)=input%delgau ; r(25)=input%tkb ; r(26)=field%efield%vslope
       r(27)=0.0 ; r(28)=0.0!r(27)=aMix_VHSE() ; r(28)=omega_VHSE()
       r(29)=input%minDistance ; r(30)=obsolete%chng ; r(31)=input%ldauMixParam ; r(32)=input%ldauSpinf
       r(33)=banddos%e_mcd_lo ; r(34)=banddos%e_mcd_up

       l(1)=input%eonly ; l(2)=input%l_useapw ; l(3)=input%secvar ; l(4)=sym%zrfs ; l(5)=input%film
       l(6)=sym%invs ; l(7)=sym%invs2 ; l(8)=input%l_bmt ; l(9)=input%l_f ; l(10)=input%cdinf
       l(11)=banddos%dos ; l(12) = hybrid%l_hybrid ; l(13)=banddos%vacdos ; l(14)=input%integ
       l(16)=input%strho ; l(17)=input%swsp ; l(18)=input%lflip 
       l(21)=input%pallst ; l(22)=sliceplot%slice ; l(23)=noco%l_soc ; l(24)=vacuum%starcoeff
       l(25)=noco%l_noco ; l(26)=noco%l_ss; l(27)=noco%l_mperp; l(28)=noco%l_constr
       l(29)=oneD%odd%d1 ; l(32)=input%ctail ; l(33)=banddos%l_orb
       l(34)=banddos%l_mcd ; l(35)=banddos%unfoldband ; l(36)=noco%l_mtNocoPot
       l(38)=field%efield%l_segmented
       l(39)=sym%symor ; l(40)=input%frcor ; l(41)=input%tria ; l(42)=field%efield%dirichlet
       l(43)=field%efield%l_dirichlet_coeff ; l(44)=input%l_coreSpec ; l(45)=input%ldauLinMix
       l(46)=noco%l_spav
    ENDIF
    !
    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,0,mpi%mpi_comm,ierr)
    hybrid%bands1=i(35) ;  input%imix=i(37);input%maxiter=i(36)
    input%gw=i(31) ; input%gw_neigd=i(32) ; hybrid%ewaldlambda=i(33) ; hybrid%lexp=i(34)
    vacuum%nstars=i(27) ; vacuum%nstm=i(28) ; oneD%odd%nq2=i(29) ; oneD%odd%nop=i(30)
    atoms%n_u=i(22) ; sym%nop2=i(23) ; sym%nsymt = i(24) 
    sliceplot%nnne=i(17) ; banddos%ndir=i(18) ; stars%mx1=i(19) ; stars%mx2=i(20) ; stars%mx3=i(21)
    input%jspins=i(12) ; vacuum%nvac=i(13) ; input%itmax=i(14) ; sliceplot%kk=i(15) ; vacuum%layers=i(16)
    stars%ng2=i(7) ; stars%ng3=i(8) ; vacuum%nmz=i(9) ; vacuum%nmzxy=i(10) ; obsolete%lepr=i(11)
     atoms%ntype=i(3) ; banddos%orbCompAtom=i(38);banddos%s_cell_x=i(40);banddos%s_cell_y=i(41);banddos%s_cell_z=i(42) ;sliceplot%iplot=i(43)
     input%coretail_lmax=i(2) ; input%kcrel=i(39)
     stars%kimax=i(25);stars%kimax2=i(26)
    !
    CALL MPI_BCAST(r,SIZE(r),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    input%minDistance=r(29) ; obsolete%chng=r(30)
    input%delgau=r(24) ; input%tkb=r(25) ; field%efield%vslope=r(26)
    cell%volint=r(19) ; hybrid%gcutm1=r(20) ; hybrid%tolerance1=r(21) 
    input%sigma=r(15) ; field%efield%zsigma=r(16); noco%mix_b=r(17); cell%vol=r(18);
    vacuum%locx(1)=r(11); vacuum%locx(2)=r(12); vacuum%locy(1)=r(13); vacuum%locy(2)=r(14)
    sliceplot%e1s=r(6) ; sliceplot%e2s=r(7) ; noco%theta=r(8); noco%phi=r(9); vacuum%tworkf=r(10)
    cell%omtil=r(1) ; cell%area=r(2) ; vacuum%delz=r(3) ; cell%z1=r(4) ; input%alpha=r(5)
    input%ldauMixParam=r(31) ; input%ldauSpinf=r(32) ; banddos%e_mcd_lo=r(33) ; banddos%e_mcd_up=r(34)
    !
    CALL MPI_BCAST(l,SIZE(l),MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    field%efield%l_dirichlet_coeff = l(43) ; input%l_useapw=l(2)
    sym%symor=l(39) ; input%frcor=l(40) ; input%tria=l(41) ; field%efield%dirichlet = l(42)
     field%efield%l_segmented=l(38)
    oneD%odd%d1=l(29) ; input%ctail=l(32)
    noco%l_noco=l(25) ; noco%l_ss=l(26) ; noco%l_mperp=l(27) ; noco%l_constr=l(28)
    input%pallst=l(21) ; sliceplot%slice=l(22) ; noco%l_soc=l(23) ; vacuum%starcoeff=l(24)
    input%strho=l(16) ; input%swsp=l(17) ; input%lflip=l(18)  
    banddos%dos=l(11) ; hybrid%l_hybrid=l(12) ; banddos%vacdos=l(13) ; banddos%l_orb=l(33) ; banddos%l_mcd=l(34)
    input%integ=l(14) 
    sym%invs=l(6) ; sym%invs2=l(7) ; input%l_bmt=l(8) ; input%l_f=l(9) ; input%cdinf=l(10)
    input%eonly=l(1)  ; input%secvar=l(3) ; sym%zrfs=l(4) ; input%film=l(5)
    field%efield%l_segmented = l(38) ; sym%symor=l(39); field%efield%dirichlet = l(40)
    field%efield%l_dirichlet_coeff = l(41) ; input%l_coreSpec=l(44) ; input%ldauLinMix=l(45)
    banddos%unfoldband=l(35)
    noco%l_mtNocoPot=l(36)
    noco%l_spav=l(46)
    !
    ! -> Broadcast the arrays:
    IF (field%efield%l_segmented) THEN
       IF (.NOT. ALLOCATED (field%efield%rhoEF))&
            &    ALLOCATE (field%efield%rhoEF(3*stars%mx1*3*stars%mx2-1,vacuum%nvac))
       n = (3*stars%mx1*3*stars%mx2-1)*vacuum%nvac
       CALL MPI_BCAST (field%efield%rhoEF,n,MPI_REAL,0,mpi%mpi_comm,ierr)
    END IF
    IF (field%efield%l_dirichlet_coeff) THEN
       IF (.NOT. ALLOCATED (field%efield%C1)) THEN
          ALLOCATE (field%efield%C1(stars%ng2-1))
          ALLOCATE (field%efield%C2(stars%ng2-1))
       END IF
       n = stars%ng2-1
       CALL MPI_BCAST (field%efield%C1,n,MPI_REAL,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST (field%efield%C2,n,MPI_REAL,0,mpi%mpi_comm,ierr)
    END IF
   
   
    n = sphhar%memd*(sphhar%nlhd+1)*sphhar%ntypsd
    CALL MPI_BCAST(sphhar%clnu,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sphhar%mlh,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sphhar%nlh,sphhar%ntypsd,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n = (sphhar%nlhd+1)*sphhar%ntypsd
    CALL MPI_BCAST(sphhar%nmem,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sphhar%llh,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%jri,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ncv,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ntypsy,atoms%nat,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%neq,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lnonsph,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lmax,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%invsat,atoms%nat,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invsatnr,atoms%nat,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ngopr,atoms%nat,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%mrot,9*sym%nop,MPI_INTEGER,0,mpi%mpi_comm,ierr)
   
    n = (2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)
    CALL MPI_BCAST(input%ellow,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%elup,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%rkmax,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%rmt,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%volmts,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%dx,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(enpara%evac,2*input%jspins*1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(enpara%evac0,2*input%jspins*1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(cell%amat,9,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(cell%bmat,9,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(cell%bbmat,9,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%taual,3*atoms%nat,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%pos,3*atoms%nat,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%tau,3*sym%nop,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    n = (atoms%lmaxd+1)*atoms%ntype*input%jspins*1
    CALL MPI_BCAST(enpara%el0,n,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    n = atoms%nlod*atoms%ntype*input%jspins
    CALL MPI_BCAST(enpara%ello0,n,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(enpara%qn_el,SIZE(enpara%qn_el),MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(enpara%qn_ello,SIZE(enpara%qn_ello),MPI_INTEGER,0,mpi%mpi_comm,ierr)


    CALL MPI_BCAST(atoms%rmsh,atoms%jmtd*atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    !
    CALL MPI_BCAST(kpts%nkpt,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%nkptf,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%bk,3*kpts%nkpt,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%wtkpt,kpts%nkpt,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%ntetra,kpts%ntet,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(kpts%voltet,kpts%ntet,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    !

    n = atoms%nat*sym%nop
    CALL MPI_BCAST(sym%invarop,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%multab,sym%nop**2,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invarind,atoms%nat,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invtab,sym%nop,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(sym%invsatnr,atoms%nat,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(vacuum%izlay,vacuum%layerd*2,MPI_INTEGER,0,mpi%mpi_comm,ierr)
 
    CALL MPI_BCAST(atoms%zatom,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(field%efield%sig_b,2,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%zelec,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ncst,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%nlo,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n =  atoms%nlod*atoms%ntype
    CALL MPI_BCAST(atoms%llo,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ulo_der,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%l_dulo,n,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    n = (atoms%llod+1)*atoms%ntype
    CALL MPI_BCAST(atoms%lo1l,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%nlol,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    n = input%jspins*atoms%ntype
    CALL MPI_BCAST(enpara%skiplo,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%alphInit,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%alph,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%beta,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%b_con,atoms%ntype*2,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%l_relax,atoms%ntype,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%l_geo,atoms%ntype,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(noco%qss,3,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%l,atoms%n_u,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%u,atoms%n_u,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%j,atoms%n_u,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%l_amf,atoms%n_u,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%atomType,atoms%n_u,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%phi,atoms%n_u,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lda_u(:)%theta,atoms%n_u,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%lapw_l,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
 
    n = 7*7*3*sym%nop
    CALL MPI_BCAST(sym%d_wgn,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(oneD%nstr1,oneD%odd%n2d,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(oneD%tau1,3*oneD%odd%nop,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    IF (oneD%odd%d1) THEN
       CALL MPI_BCAST(oneD%tau1,3*oneD%odd%nop,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%mrot1,9*oneD%odd%nop,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%kv1,2*oneD%odd%n2d,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%ngopr1,atoms%nat,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%igfft1,oneD%odd%nn2d*2,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%pgfft1,oneD%odd%nn2d,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
       n = (2*stars%mx3 + 1)*(2*ONED%ODD%M +1)
       CALL MPI_BCAST(oneD%ig1,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%invtab1,oneD%odd%nop,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(oneD%multab1,2*oneD%odd%nop,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    ENDIF
    !--- HF<
    CALL MPI_BCAST(kpts%nkpt3,3,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    IF(hybrid%l_hybrid) THEN
       CALL MPI_BCAST(hybrid%lcutwf,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(hybrid%select1,4*atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(hybrid%lcutm1,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    END IF
    !--- HF>

    IF (mpi%irank>0) THEN
       ALLOCATE(noco%socscale(atoms%ntype))
    ENDIF
    CALL MPI_BCAST(noco%socscale,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)

    IF(input%l_inpXML) THEN
       n = dimension%nstd*atoms%ntype
       CALL MPI_BCAST(atoms%numStatesProvided,atoms%ntype,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(atoms%coreStateOccs,2*n,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(atoms%coreStateNprnc,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(atoms%coreStateKappa,n,MPI_INTEGER,0,mpi%mpi_comm,ierr)

       CALL MPI_BCAST(kpts%posScale,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(kpts%numSpecialPoints,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(kpts%specialPoints,3*kpts%numSpecialPoints,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    END IF
 
    IF(banddos%unfoldband) THEN
       IF(mpi%irank.NE.0) THEN
          ALLOCATE(kpts%sc_list(13,kpts%nkpt))
       END IF
       CALL MPI_BCAST(kpts%sc_list,13*kpts%nkpt,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    END IF

    RETURN
#endif
  END SUBROUTINE mpi_bc_all

  SUBROUTINE priv_mpi_bc_stars(mpi,stars)
    USE m_types
    USE m_mpi_bc_tool
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(in)::mpi
    TYPE(t_stars),INTENT(inout)::stars
    include 'mpif.h'

    INTEGER :: i(18),ierr,ft2_gf_dim
    
    i(1)=stars%ng2 ; i(2)=stars%ng3 ; i(3)=stars%mx1 ; i(4)=stars%mx2 ; i(5)=stars%mx3
    i(6) = stars%kimax ; i(7) = stars%kimax2 ; i(8)=size(stars%ft2_gfy)
    i(9) = stars%kq1_fft;i(10) = stars%kq2_fft;i(11) = stars%kq3_fft;i(12)=stars%kmxq_fft
    i(13)= stars%kxc1_fft;i(14)= stars%kxc2_fft;i(15)= stars%kxc3_fft
    i(16)= stars%ng3_fft;i(17)= stars%kmxxc_fft;i(18)= stars%nxc3_fft
    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,0,mpi%mpi_comm,ierr)
    stars%ng2=i(1) ; stars%ng3=i(2) ; stars%mx1=i(3) ; stars%mx2=i(4) ; stars%mx3=i(5)
    stars%kimax=i(6) ;  stars%kimax2=i(7) ; ft2_gf_dim=i(8)
    stars%kq1_fft=i(9);stars%kq2_fft=i(10); stars%kq3_fft=i(11);stars%kmxq_fft=i(12)
    stars%kxc1_fft=i(13);stars%kxc2_fft=i(14);stars%kxc3_fft=i(15)
    stars%ng3_fft=i(16);stars%kmxxc_fft=i(17);stars%nxc3_fft=i(18)


    CALL mpi_bc(stars%ustep,0,mpi%mpi_comm)
    CALL mpi_bc(stars%ig2,0,mpi%mpi_comm)
    CALL mpi_bc(stars%ig,0,mpi%mpi_comm)
    CALL mpi_bc(stars%rgphs,0,mpi%mpi_comm)
    CALL mpi_bc(stars%sk3,0,mpi%mpi_comm)
    CALL mpi_bc(stars%kv2,0,mpi%mpi_comm)
    CALL mpi_bc(stars%kv3,0,mpi%mpi_comm)
    CALL mpi_bc(stars%nstr,0,mpi%mpi_comm)
    CALL mpi_bc(stars%nstr2,0,mpi%mpi_comm)
    CALL mpi_bc(stars%igfft,0,mpi%mpi_comm)
    CALL mpi_bc(stars%pgfft,0,mpi%mpi_comm)
    CALL mpi_bc(stars%igfft2,0,mpi%mpi_comm)
    CALL mpi_bc(stars%sk2,0,mpi%mpi_comm)
    CALL mpi_bc(stars%phi2,0,mpi%mpi_comm)
    CALL mpi_bc(stars%ufft,0,mpi%mpi_comm)
    CALL mpi_bc(stars%igq2_fft,0,mpi%mpi_comm)
    CALL mpi_bc(stars%igq_fft,0,mpi%mpi_comm)
    CALL mpi_bc(stars%pgfft2,0,mpi%mpi_comm)
    CALL mpi_bc(stars%ft2_gfx,0,mpi%mpi_comm)
    CALL mpi_bc(stars%ft2_gfy,0,mpi%mpi_comm)
  END SUBROUTINE priv_mpi_bc_stars

END MODULE m_mpi_bc_all
