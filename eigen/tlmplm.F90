MODULE m_tlmplm
  use m_judft
  IMPLICIT NONE
  !*********************************************************************
  !     sets up the local Hamiltonian, i.e. the Hamiltonian in the
  !     l',m',l,m,u- basis which is independent from k!
  !*********************************************************************
CONTAINS
  SUBROUTINE tlmplm(n,sphhar,atoms,enpara,&
       jspin,jsp,mpi,v,input,td,ud)
    USE m_constants
    USE m_intgr, ONLY : intgr3
    USE m_genMTBasis
    USE m_tlo
    USE m_gaunt, ONLY: gaunt1
    USE m_types
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_potden),INTENT(IN)    :: v
    TYPE(t_tlmplm),INTENT(INOUT) :: td
    TYPE(t_usdus),INTENT(INOUT)  :: ud

    INTEGER, INTENT (IN) :: n,jspin !atom index,physical spin&spin index for data

    REAL, ALLOCATABLE   :: dvd(:,:),dvu(:,:),uvd(:,:),uvu(:,:),f(:,:,:,:),g(:,:,:,:),x(:),flo(:,:,:)
    INTEGER,ALLOCATABLE :: indt(:)
    REAL,ALLOCATABLE    :: vr0(:,:)
   

    COMPLEX  :: cil
    REAL     :: temp
    INTEGER i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmpl,lmplm,lmx,lmxx,lp,info,in
    INTEGER lp1,lpl ,mem,mems,mp,mu,nh,na,m,nsym,s,i_u,jspin1,jspin2

    ALLOCATE( dvd(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd ))
    ALLOCATE( dvu(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd ))
    ALLOCATE( uvd(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd ))
    ALLOCATE( uvu(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd ))
    ALLOCATE( f(atoms%jmtd,2,0:atoms%lmaxd,2),g(atoms%jmtd,2,0:atoms%lmaxd,2),x(atoms%jmtd))

    ALLOCATE( flo(atoms%jmtd,2,atoms%nlod))
    ALLOCATE( indt(0:SIZE(td%tuu,1)-1))
    ALLOCATE( vr0(SIZE(v%mt,1),0:SIZE(v%mt,2)-1))


    jsp=jspin
    vr0=v%mt(:,:,n,jsp)
    IF (jsp<3) vr0(:,0)=0.0

    DO i=MERGE(1,jspin,jspin>2),MERGE(2,jspin,jspin>2)
       CALL genMTBasis(atoms,enpara,v,mpi,n,i,ud,f(:,:,:,i),g(:,:,:,i),flo)
    ENDDO
    IF (jspin>2) THEN
       jspin1=1
       jspin2=2
    ELSE
       jspin1=jspin;jspin2=jspin
    END IF
    na=SUM(atoms%neq(:n-1))+1
    nsym = atoms%ntypsy(na)
    nh = sphhar%nlh(nsym)
    !
    !--->    generate the irreducible integrals (u(l'):v(lamda,nu:u(l))
    !--->    for l' .ge. l, but only those that will contribute
    !
    DO lp = 0,atoms%lmax(n)
       lp1 = (lp* (lp+1))/2
       DO l = 0,lp
          lpl = lp1 + l
          !--->    loop over non-spherical components of the potential: must
          !--->    satisfy the triangular conditions and that l'+l+lamda even
          !--->    (conditions from the gaunt coefficient)
          DO lh = MERGE(1,0,jspin<3), nh
             lamda = sphhar%llh(lh,nsym)
             lmin = lp - l
             lmx = lp + l
             IF ((mod(lamda+lmx,2).EQ.1) .OR. (lamda.LT.lmin) .OR. (lamda.GT.lmx)) THEN
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

    td%tuu(0:,n,jsp) = cmplx(0.0,0.0)
    td%tdd(0:,n,jsp) = cmplx(0.0,0.0)
    td%tud(0:,n,jsp) = cmplx(0.0,0.0)
    td%tdu(0:,n,jsp) = cmplx(0.0,0.0)
    indt=0
    !--->    generate the various t(l'm',lm) matrices for l'm'.ge.lm
    !--->    loop over l'm'
    DO lp = 0,atoms%lmax(n)
       lp1 = (lp* (lp+1))/2
       DO mp = -lp,lp
          lmp = lp* (lp+1) + mp
          lmpl = (lmp* (lmp+1))/2
          !--->    loop over lattice harmonics
          DO lh = MERGE(1,0,jspin<3), nh
             lamda = sphhar%llh(lh,nsym)
             lmin0 = abs(lp-lamda)
             IF (lmin0.GT.lp) CYCLE
             !-->     ensure l+l'+lamda even
             lmxx = lp - mod(lamda,2)
             mems = sphhar%nmem(lh,nsym)
             DO mem = 1,mems
                mu = sphhar%mlh(mem,lh,nsym)
                m = mp - mu
                lmin = max(lmin0,abs(m))
                l2 = abs(lmxx-lmin)
                lmin = lmin + mod(l2,2)
                DO l = lmin,lmxx,2
                   lm = l* (l+1) + m
                   IF (lm.GT.lmp) CYCLE
                   lpl = lp1 + l
                   lmplm = lmpl + lm
                   cil = ((ImagUnit** (l-lp))*sphhar%clnu(mem,lh,nsym))*&
                        gaunt1(lp,lamda,l,mp,mu,m,atoms%lmaxd)
                   td%tuu(lmplm,n,jsp) = td%tuu(lmplm,n,jsp) + cil*uvu(lpl,lh)
                   td%tdd(lmplm,n,jsp) = td%tdd(lmplm,n,jsp) + cil*dvd(lpl,lh)
                   td%tud(lmplm,n,jsp) = td%tud(lmplm,n,jsp) + cil*uvd(lpl,lh)
                   td%tdu(lmplm,n,jsp) = td%tdu(lmplm,n,jsp) + cil*dvu(lpl,lh)
                   indt(lmplm) = 1
                END DO
             END DO
          END DO
       END DO
    END DO
    !--->    set up mapping array
    DO lp = 0,atoms%lmax(n)
       DO mp = -lp,lp
          lmp = lp* (lp+1) + mp
          DO l = 0,atoms%lmax(n)
             DO m = -l,l
                lm = l* (l+1) + m
                IF (lmp.GE.lm) THEN
                   lmplm = (lmp* (lmp+1))/2 + lm
                   IF (indt(lmplm).NE.0) THEN
                      td%ind(lmp,lm,n,jsp) = lmplm
                   ELSE
                      td%ind(lmp,lm,n,jsp) = -9999
                   END IF
                ELSE
                   lmplm = (lm* (lm+1))/2 + lmp
                   IF (indt(lmplm).NE.0) THEN
                      td%ind(lmp,lm,n,jsp) = -lmplm
                   ELSE
                      td%ind(lmp,lm,n,jsp) = -9999
                   END IF
                END IF
             END DO
          END DO
       END DO
    ENDDO

    !
    !--->   set up the t-matrices for the local orbitals,
    !--->   if there are any
    IF (atoms%nlo(n).GE.1) THEN
       IF (jspin>3) call judft_error("l_mtnocoPot=T and LOs not supported.")
          CALL tlo(atoms,sphhar,jspin,jsp,n,enpara,1,input,v%mt(1,0,n,jsp),&
               na,flo,f(:,:,:,jspin),g(:,:,:,jspin),ud, ud%uuilon(:,:,jspin),ud%duilon(:,:,jspin),ud%ulouilopn(:,:,:,jspin), td)

    ENDIF
  END SUBROUTINE tlmplm
END MODULE m_tlmplm
