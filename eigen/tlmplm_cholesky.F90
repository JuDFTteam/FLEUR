MODULE m_tlmplm_cholesky

  IMPLICIT NONE
  !*********************************************************************
  !     sets up the local Hamiltonian, i.e. the Hamiltonian in the
  !     l',m',l,m,u- basis which is independent from k!
  !     shifts this local Hamiltonian to make it positive definite
  !     and does a cholesky decomposition
  !*********************************************************************
  CONTAINS
    SUBROUTINE tlmplm_cholesky(sphhar,atoms,noco,enpara,&
         jspin,jsp,mpi,v,input,td,ud)

      USE m_intgr, ONLY : intgr3
      USE m_genMTBasis
      USE m_tlo
      USE m_gaunt, ONLY: gaunt1,gaunt2
      USE m_types
      USE m_radovlp
      IMPLICIT NONE
      TYPE(t_mpi),INTENT(IN)      :: mpi
      TYPE(t_noco),INTENT(IN)     :: noco
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_enpara),INTENT(IN)   :: enpara
      !     ..
      !     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspin,jsp !physical spin&spin index for data
      !     ..
      TYPE(t_potden),INTENT(IN)   :: v
      TYPE(t_tlmplm),INTENT(INOUT) :: td
      TYPE(t_usdus),INTENT(INOUT)  :: ud
      
      !     ..
      !     .. Local Scalars ..
      COMPLEX cil
      COMPLEX,PARAMETER::ci=cmplx(0.,1.)
      REAL temp
      INTEGER i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmpl,lmplm,lmx,lmxx,lp,info,in
      INTEGER lp1,lpl ,mem,mems,mp,mu,n,nh,na,m,nsym,s,i_u
      LOGICAL l_write,OK
      !     ..
      !     .. Local Arrays ..
      REAL vr0(size(v%mt,1),0:size(v%mt,2)-1,size(v%mt,3))
      REAL dvd(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL dvu(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL uvd(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL uvu(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL f(atoms%jmtd,2,0:atoms%lmaxd),g(atoms%jmtd,2,0:atoms%lmaxd),x(atoms%jmtd)
      REAL flo(atoms%jmtd,2,atoms%nlod)
      INTEGER:: indt(0:SIZE(td%tuu,1)-1)

      !for constraint
      REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)
      COMPLEX :: c
      
      REAL,PARAMETER:: e_shift_min=0.5
      REAL,PARAMETER:: e_shift_max=65.0
    
    vr0=v%mt(:,:,:,jsp)
    vr0(:,0,:)=0.0
    !     ..e_shift
    td%e_shift(:,jsp)=e_shift_min

   
     IF (noco%l_constr) THEN
       ALLOCATE(uun21(0:atoms%lmaxd,atoms%ntype),udn21(0:atoms%lmaxd,atoms%ntype),&
            dun21(0:atoms%lmaxd,atoms%ntype),ddn21(0:atoms%lmaxd,atoms%ntype) )
       CALL rad_ovlp(atoms,ud,input,v%mt,enpara%el0, uun21,udn21,dun21,ddn21)
    ENDIF
    
    l_write = mpi%irank==0

    td%tdulo(:,:,:,jsp) = cmplx(0.0,0.0)
    td%tuulo(:,:,:,jsp) = cmplx(0.0,0.0)
    td%tuloulo(:,:,:,jsp) = cmplx(0.0,0.0)

    td%h_off=0.0
!$    l_write=.false.
!$    call gaunt2(atoms%lmaxd)
!$OMP PARALLEL DO DEFAULT(NONE)&
!$OMP PRIVATE(indt,dvd,dvu,uvd,uvu,f,g,x,flo)&
!$OMP PRIVATE(cil,temp,i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmpl)&
!$OMP PRIVATE(lmplm,lmx,lmxx,lp,lp1,lpl,m,mem,mems,mp,mu,n,nh)&
!$OMP PRIVATE(nsym,na,OK,s,in,info,c)&
!$OMP SHARED(atoms,jspin,jsp,sphhar,enpara,td,ud,l_write,v,mpi,input,vr0)&
!$OMP SHARED(noco,uun21,udn21,dun21,ddn21)
    DO  n = 1,atoms%ntype
       na=sum(atoms%neq(:n-1))+1
       OK=.FALSE.

       cholesky_loop:DO WHILE(.NOT.OK)
          td%h_loc(:,:,n,jsp)=0.0
          OK=.TRUE.
          !
          !--->    generate the wavefunctions for each l
          !
          CALL genMTBasis(atoms,enpara,v,mpi,n,jspin,l_write,ud,f,g,flo)
          
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
                DO lh = 1, nh
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
                         x(i) = (f(i,1,lp)*f(i,1,l)+f(i,2,lp)*f(i,2,l))* vr0(i,lh,n)
                      END DO
                      CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                      uvu(lpl,lh) = temp
                      DO i = 1,atoms%jri(n)
                         x(i) = (g(i,1,lp)*f(i,1,l)+g(i,2,lp)*f(i,2,l))* vr0(i,lh,n)
                      END DO
                      CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                      dvu(lpl,lh) = temp
                      DO i = 1,atoms%jri(n)
                         x(i) = (f(i,1,lp)*g(i,1,l)+f(i,2,lp)*g(i,2,l))* vr0(i,lh,n)
                      END DO
                      CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                      uvd(lpl,lh) = temp
                      DO i = 1,atoms%jri(n)
                         x(i) = (g(i,1,lp)*g(i,1,l)+g(i,2,lp)*g(i,2,l))* vr0(i,lh,n)
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
                DO lh = 1, nh
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
                         cil = ((ci** (l-lp))*sphhar%clnu(mem,lh,nsym))*&
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

          s=atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+1
          !Setup local hamiltonian
          DO lmp=0,atoms%lnonsph(n)*(atoms%lnonsph(n)+2)
             lp=FLOOR(SQRT(1.0*lmp))
             mp=lmp-lp*(lp+1)
             IF (lp>atoms%lmax(n).OR.ABS(mp)>lp) STOP "BUG"
             !--->             loop over l,m
             DO l = 0,atoms%lnonsph(n)
                DO m = -l,l
                   lm = l* (l+1) + m
                   in = td%ind(lmp,lm,n,jsp)
                   IF (in/=-9999) THEN
                      IF (in>=0) THEN
                         td%h_loc(lm,lmp,n,jsp)    = CONJG(td%tuu(in,n,jsp))
                         td%h_loc(lm+s,lmp,n,jsp)  = CONJG(td%tud(in,n,jsp))
                         td%h_loc(lm,lmp+s,n,jsp)  = CONJG(td%tdu(in,n,jsp))
                         td%h_loc(lm+s,lmp+s,n,jsp)= CONJG(td%tdd(in,n,jsp))
                      ELSE
                         td%h_loc(lm,lmp,n,jsp)    = td%tuu(-in,n,jsp)
                         td%h_loc(lm+s,lmp,n,jsp)  = td%tdu(-in,n,jsp)
                         td%h_loc(lm,lmp+s,n,jsp)  = td%tud(-in,n,jsp)
                         td%h_loc(lm+s,lmp+s,n,jsp)= td%tdd(-in,n,jsp)
                      END IF
                   END IF
                END DO
             END DO
          ENDDO

          
          !Include contribution from LDA+U
          DO i_u=1,atoms%n_u
             IF (n.NE.atoms%lda_u(i_u)%atomtype) CYCLE
             !Found a "U" for this atom type
             l=atoms%lda_u(i_u)%l
             lp=atoms%lda_u(i_u)%l
             DO m = -l,l
                lm = l* (l+1) + m
                DO mp = -lp,lp
                   lmp = lp* (lp+1) + mp
                   td%h_loc(lm,lmp,n,jsp)     =td%h_loc(lm,lmp,n,jsp) + v%mmpMat(m,mp,i_u,jsp)
                   td%h_loc(lm+s,lmp+s,n,jsp) =td%h_loc(lm+s,lmp+s,n,jsp)+ v%mmpMat(m,mp,i_u,jsp)*ud%ddn(lp,n,jsp)
                ENDDO
             ENDDO
          END DO
          !Now add diagonal contribution to matrices
          DO l = 0,atoms%lmax(n)
             DO  m = -l,l
                lm = l* (l+1) + m
                lmplm = (lm* (lm+3))/2
                td%tuu(lmplm,n,jsp)=td%tuu(lmplm,n,jsp) + enpara%el0(l,n,jsp)
                td%tdd(lmplm,n,jsp)=td%tdd(lmplm,n,jsp) + enpara%el0(l,n,jsp)*ud%ddn(l,n,jsp)
                td%tud(lmplm,n,jsp)=td%tud(lmplm,n,jsp) + 0.5
                td%tdu(lmplm,n,jsp)=td%tdu(lmplm,n,jsp) + 0.5
             ENDDO
          ENDDO
          
          !Create Cholesky decomposition of local hamiltonian
          
          !--->    Add diagonal terms to make matrix positive definite
          DO lp = 0,atoms%lnonsph(n)
             DO mp = -lp,lp
                lmp = lp* (lp+1) + mp
                td%h_loc(lmp,lmp,n,jsp)=td%e_shift(n,jsp)+td%h_loc(lmp,lmp,n,jsp)
                td%h_loc(lmp+s,lmp+s,n,jsp)=td%e_shift(n,jsp)*ud%ddn(lp,n,jsp)+td%h_loc(lmp+s,lmp+s,n,jsp)
             END DO
          END DO
          IF (lmp+1.ne.s) call judft_error("BUG in tlmpln_cholesky")
          !Perform cholesky decomposition
          info=0
          CALL zpotrf("L",2*s,td%h_loc(:,:,n,jsp),SIZE(td%h_loc,1),info)

          !Upper part to zero
          DO l=0,2*s-1
             DO lp=0,l-1
                td%h_loc(lp,l,n,jsp)=0.0
             ENDDO
          ENDDO
         
          IF (info.NE.0) THEN
             td%e_shift(n,jsp)=td%e_shift(n,jsp)*2.0
             PRINT *,"Potential shift to small, increasing the value to:",td%e_shift(n,jsp)
             IF (td%e_shift(n,jsp)>e_shift_max) THEN
                 CALL judft_error("Potential shift at maximum")
             ENDIF
             OK=.FALSE.
             CYCLE cholesky_loop
          ENDIF

          !
          !--->   set up the t-matrices for the local orbitals,
          !--->   if there are any
          IF (atoms%nlo(n).GE.1) THEN
             CALL tlo(atoms,sphhar,jspin,jsp,n,enpara,1,input,v%mt(1,0,n,jsp),&
                  na,flo,f,g,ud, ud%uuilon(:,:,jspin),ud%duilon(:,:,jspin),ud%ulouilopn(:,:,:,jspin), td)
             
          ENDIF

          !If we do  a constraint calculation, we have to calculate the
          !local spin off-diagonal contributions
          s=atoms%lnonsph(n)+1
          !first ispin=2,jspin=1 case
          IF (noco%l_constr) THEN
             DO l=0,atoms%lnonsph(n)
                c=(-0.5)*CMPLX(noco%b_con(1,n),noco%b_con(2,n))
                td%h_off(l  ,l  ,n,1)     =td%h_off(l  ,l  ,n,1) + uun21(l,n)*c
                td%h_loc(l  ,l+s,n,1)     =td%h_off(l  ,l+s,n,1) + udn21(l,n)*c
                td%h_loc(l+s,l  ,n,1)     =td%h_off(l+s,l  ,n,1) + dun21(l,n)*c
                td%h_loc(l+s,l+s,n,1)     =td%h_off(l+s,l+s,n,1) + ddn21(l,n)*c
             ENDDO
          ENDIF
          
          
          !then ispin=2,jspin=1 case
          IF (noco%l_constr) THEN
             DO l=0,atoms%lnonsph(n)
                c=(-0.5)*CMPLX(noco%b_con(1,n),-noco%b_con(2,n))
                td%h_off(l  ,l  ,n,2)     =td%h_off(l  ,l  ,n,2) + uun21(l,n)*c
                td%h_loc(l  ,l+s,n,2)     =td%h_off(l  ,l+s,n,2) + udn21(l,n)*c
                td%h_loc(l+s,l  ,n,2)     =td%h_off(l+s,l  ,n,2) + dun21(l,n)*c
                td%h_loc(l+s,l+s,n,2)     =td%h_off(l+s,l+s,n,2) + ddn21(l,n)*c
             ENDDO
          ENDIF
          
          
       ENDDO cholesky_loop
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE tlmplm_cholesky


   
  
  
END MODULE m_tlmplm_cholesky
