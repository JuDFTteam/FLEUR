MODULE m_tlmplm
  USE m_judft
  IMPLICIT NONE
  !*********************************************************************
  !     sets up the local Hamiltonian, i.e. the Hamiltonian in the
  !     l',m',l,m,u- basis which is independent from k!
  !*********************************************************************
CONTAINS
  SUBROUTINE tlmplm(n,sphhar,atoms,sym,enpara,nococonv,&
       jspin1,jspin2,jsp,fmpi,v,vx,input,hub1inp,hub1data,td,ud,alpha_hybrid,lh0,one,v1)
    USE m_constants
    USE m_intgr, ONLY : intgr3
    USE m_genMTBasis
    USE m_tlo
    USE m_gaunt, ONLY: gaunt1
    USE m_types
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_enpara),INTENT(IN)    :: enpara
    TYPE(t_nococonv),INTENT(IN)  :: nococonv
    TYPE(t_mpi),INTENT(IN)       :: fmpi
    TYPE(t_potden),INTENT(IN)    :: v,vx
    TYPE(t_hub1inp),INTENT(IN)   :: hub1inp
    TYPE(t_hub1data),INTENT(IN)  :: hub1data
    TYPE(t_tlmplm),INTENT(INOUT) :: td
    TYPE(t_usdus),INTENT(INOUT)  :: ud

    INTEGER, INTENT (IN) :: n,jspin1,jspin2,jsp,lh0 !atom index,physical spin&spin index for data
    REAL,INTENT(IN)      :: alpha_hybrid

    COMPLEX, INTENT(IN) :: one

    TYPE(t_potden), OPTIONAL, INTENT(IN) :: v1

    REAL, ALLOCATABLE   :: dvd(:,:),dvu(:,:),uvd(:,:),uvu(:,:),f(:,:,:,:),g(:,:,:,:),x(:),flo(:,:,:,:)
    REAL,ALLOCATABLE    :: vr0(:,:)

    COMPLEX  :: cil
    REAL     :: temp
    INTEGER i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmpl,lmplm,lmx,lmxx,lp,info,in
    INTEGER lp1,lpl ,mem,mems,mp,mu,nh,na,m,nsym,s,i_u,lplmax
    LOGICAL l_remove

    lplmax = atoms%lmaxd*(atoms%lmaxd+3)/2

    ALLOCATE( dvd(0:lplmax,0:sphhar%nlhd ));dvd=0.0
    ALLOCATE( dvu(0:lplmax,0:sphhar%nlhd ));dvu=0.0
    ALLOCATE( uvd(0:lplmax,0:sphhar%nlhd ));uvd=0.0
    ALLOCATE( uvu(0:lplmax,0:sphhar%nlhd ));uvu=0.0
    ALLOCATE( f(atoms%jmtd,2,0:atoms%lmaxd,2),g(atoms%jmtd,2,0:atoms%lmaxd,2),x(atoms%jmtd))

    ALLOCATE( flo(atoms%jmtd,2,atoms%nlod,2))
    ALLOCATE( vr0(SIZE(v%mt,1),0:SIZE(v%mt,2)-1))

    IF (.NOT.PRESENT(v1)) THEN
        vr0=v%mt(:,:,n,jsp)
        IF (jsp<3) THEN
            vr0(:,0)=0.0
            if (alpha_hybrid.ne.0) vr0=vr0-alpha_hybrid*vx%mt(:,:,n,jsp)
        ELSE
            ! TODO: Shouldn't we correct v0 for the l=0-factor?
            vr0(:,0)=vr0(:,0)-0.5*nococonv%b_con(jsp-2,n) !Add constraining field
        ENDIF
    ELSE
        vr0=v1%mt(:, :, n, jsp)
    END IF

    DO i=MIN(jspin1,jspin2),MAX(jspin1,jspin2)
       CALL genMTBasis(atoms,enpara,v,fmpi,n,i,ud,f(:,:,:,i),g(:,:,:,i),flo(:,:,:,i),hub1data=hub1data)
    ENDDO
    na=SUM(atoms%neq(:n-1))+1
    nsym = sym%ntypsy(na)
    nh = sphhar%nlh(nsym)
    !
    !--->    generate the irreducible integrals (u(l'):v(lamda,nu:u(l))
    !--->    for l' .ge. l, but only those that will contribute
    !
    DO lp = 0,atoms%lmax(n)
       lp1 = (lp*(lp+1))/2
       DO l = 0, lp
          lpl = lp1 + l
          !----------------------------------------------------------------------------
          ! Remove non-spherical components for the orbitals treated with DFT+Hubbard-1
          !----------------------------------------------------------------------------
          l_remove=.FALSE.
          IF(l.EQ.lp.AND.hub1inp%l_nonsphDC) THEN
             DO i = atoms%n_u+1, atoms%n_u+atoms%n_hia
                IF(atoms%lda_u(i)%atomType.EQ.n.AND.atoms%lda_u(i)%l.EQ.l) l_remove=.TRUE.
             ENDDO
          ENDIF
          !--->    loop over non-spherical components of the potential: must
          !--->    satisfy the triangular conditions and that l'+l+lamda even
          !--->    (conditions from the gaunt coefficient)
          DO lh = lh0, nh
             lamda = sphhar%llh(lh,nsym)
             lmin = lp - l
             lmx = lp + l
             IF ((MOD(lamda+lmx,2).EQ.1) .OR. (lamda.LT.lmin) .OR. (lamda.GT.lmx) .OR. l_remove) THEN
                uvu(lpl,lh) = 0.0
                dvd(lpl,lh) = 0.0
                uvd(lpl,lh) = 0.0
                dvu(lpl,lh) = 0.0
             ELSE
                DO i = 1,atoms%jri(n)
                   x(i) = (f(i,1,lp,jspin1)*f(i,1,l,jspin2)+f(i,2,lp,jspin1)*f(i,2,l,jspin2))* vr0(i,lh)
                END DO
                CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                uvu(lpl,lh) = temp
                DO i = 1,atoms%jri(n)
                   x(i) = (g(i,1,lp,jspin1)*f(i,1,l,jspin2)+g(i,2,lp,jspin1)*f(i,2,l,jspin2))* vr0(i,lh)
                END DO
                CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                dvu(lpl,lh) = temp
                DO i = 1,atoms%jri(n)
                   x(i) = (f(i,1,lp,jspin1)*g(i,1,l,jspin2)+f(i,2,lp,jspin1)*g(i,2,l,jspin2))* vr0(i,lh)
                END DO
                CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                uvd(lpl,lh) = temp
                DO i = 1,atoms%jri(n)
                   x(i) = (g(i,1,lp,jspin1)*g(i,1,l,jspin2)+g(i,2,lp,jspin1)*g(i,2,l,jspin2))* vr0(i,lh)
                END DO
                CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                dvd(lpl,lh) = temp
             END IF
          END DO
       END DO
    END DO


    !--->    generate the various t(l'm',lm) matrices for l'm'.ge.lm
    !--->    loop over l'm'
    s=td%h_loc2(n)
    DO lp = 0,atoms%lmax(n)
       lp1 = (lp*(lp+1))/2
       DO mp = -lp,lp
          lmp = lp* (lp+1) + mp
          !lmpl = (lmp* (lmp+1))/2
          !--->    loop over lattice harmonics
          DO lh = lh0, nh
             lamda = sphhar%llh(lh,nsym)
             lmin0 = ABS(lp-lamda)
             IF (lmin0.GT.lp) CYCLE
             !-->     ensure l+l'+lamda even
             lmxx = lp - MOD(lamda,2)
             mems = sphhar%nmem(lh,nsym)
             DO mem = 1,mems
                mu = sphhar%mlh(mem,lh,nsym)
                m = mp - mu
                lmin = MAX(lmin0,ABS(m))
                l2 = ABS(lmxx-lmin)
                lmin = lmin + MOD(l2,2)
                DO l = lmin,lmxx,2
                   lm = l* (l+1) + m
                   IF (lm.GT.lmp) CYCLE
                   lpl = lp1 + l
                   !lmplm = lmpl + lm
                   cil = ((ImagUnit** (l-lp))*sphhar%clnu(mem,lh,nsym))*&
                        gaunt1(lp,lamda,l,mp,mu,m,atoms%lmaxd)
                   td%h_loc(lmp,lm,n,jspin1,jspin2)    =  td%h_loc(lmp,lm,n,jspin1,jspin2)+ one*cil*uvu(lpl,lh)
                   td%h_loc(lmp+s,lm,n,jspin1,jspin2)  =  td%h_loc(lmp+s,lm,n,jspin1,jspin2)+ one*cil*dvu(lpl,lh)
                   td%h_loc(lmp,lm+s,n,jspin1,jspin2)  =  td%h_loc(lmp,lm+s,n,jspin1,jspin2)+ one*cil*uvd(lpl,lh)
                   td%h_loc(lmp+s,lm+s,n,jspin1,jspin2)    =  td%h_loc(lmp+s,lm+s,n,jspin1,jspin2)+ one*cil*dvd(lpl,lh)
                   IF (lm.NE.lmp) THEN
                      td%h_loc(lm,lmp,n,jspin1,jspin2)    =  td%h_loc(lm,lmp,n,jspin1,jspin2)+ one*CONJG(cil*uvu(lpl,lh))
                      td%h_loc(lm+s,lmp,n,jspin1,jspin2)  =  td%h_loc(lm+s,lmp,n,jspin1,jspin2)+ one*CONJG(cil*uvd(lpl,lh))
                      td%h_loc(lm,lmp+s,n,jspin1,jspin2)  =  td%h_loc(lm,lmp+s,n,jspin1,jspin2)+ one*CONJG(cil*dvu(lpl,lh))
                      td%h_loc(lm+s,lmp+s,n,jspin1,jspin2)    =  td%h_loc(lm+s,lmp+s,n,jspin1,jspin2)+ one*CONJG(cil*dvd(lpl,lh))
                   ENDIF
                END DO
             END DO
          END DO
       END DO
    END DO

    !
    !--->   set up the t-matrices for the local orbitals,
    !--->   if there are any
    IF (atoms%nlo(n).GE.1) THEN
       CALL tlo(atoms,sym,sphhar,jspin1,jspin2,jsp,n,enpara,lh0,input,vr0,&
            na,flo,f,g,ud, td, one)
    ENDIF
  END SUBROUTINE tlmplm
END MODULE m_tlmplm
