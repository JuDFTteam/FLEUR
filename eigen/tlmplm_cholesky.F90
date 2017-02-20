MODULE m_tlmplm_cholesky
  IMPLICIT NONE
  !*********************************************************************
  !     sets up the local Hamiltonian, i.e. the Hamiltonian in the
  !     l',m',l,m,u- basis which is independent from k!
  !     shifts this local Hamiltonian to make it positive definite
  !     and does a cholesky decomposition
  !*********************************************************************
  CONTAINS
    SUBROUTINE tlmplm_cholesky(sphhar,atoms,dimension,enpara,&
         jspin,jsp,mpi, vr,input, td,ud)

      USE m_intgr, ONLY : intgr3
      USE m_radflo
      USE m_radfun
      USE m_tlo
      USE m_gaunt, ONLY: gaunt1,gaunt2
      USE m_types
      IMPLICIT NONE
      TYPE(t_mpi),INTENT(IN)      :: mpi
      TYPE(t_dimension),INTENT(IN):: dimension
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_enpara),INTENT(IN)   :: enpara
      !     ..
      !     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspin,jsp !physical spin&spin index for data
      !     ..
      !     .. Array Arguments ..
      
      REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)   ! this is for the
      TYPE(t_tlmplm),INTENT(INOUT) :: td
      TYPE(t_usdus),INTENT(INOUT)  :: ud
      
      !     ..
      !     .. Local Scalars ..
      COMPLEX cil
      COMPLEX,PARAMETER::ci=cmplx(0.,1.)
      REAL temp,wronk
      INTEGER i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmpl,lmplm,lmx,lmxx,lp,info,in
      INTEGER lp1,lpl ,mem,mems,mp,mu,n,nh,noded,nodeu ,na,m,nsym,s
      LOGICAL l_write,ok
      !     ..
      !     .. Local Arrays ..
      REAL vr0(size(vr,1),0:size(vr,2)-1,size(vr,3))
      REAL dvd(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL dvu(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL uvd(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL uvu(0:atoms%lmaxd*(atoms%lmaxd+3)/2,0:sphhar%nlhd )
      REAL f(atoms%jmtd,2,0:atoms%lmaxd),g(atoms%jmtd,2,0:atoms%lmaxd),x(atoms%jmtd)
      REAL flo(atoms%jmtd,2,atoms%nlod)
      REAL uuilon(atoms%nlod,atoms%ntype),duilon(atoms%nlod,atoms%ntype)
      REAL ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
      INTEGER:: indt(0:dimension%lmplmd)

    REAL,PARAMETER:: e_shift_min=64.0
    REAL,PARAMETER:: e_shift_max=20000000.0
    
    vr0=vr
    vr0(:,0,:)=0.0
    !     ..e_shift
    td%e_shift=e_shift_min
    OK=.false.
    DO WHILE(.not.OK)
       td%h_loc=0.0
       OK=.true.
       td%tdulo(:,:,:,jsp) = cmplx(0.0,0.0)
       td%tuulo(:,:,:,jsp) = cmplx(0.0,0.0)
       td%tuloulo(:,:,:,jsp) = cmplx(0.0,0.0)

       !
       !--->    generate the wavefunctions for each l
       !
       l_write=mpi%irank==0
!!$    l_write=.false.
!!$    call gaunt2(atoms%lmaxd)
!!$OMP PARALLEL DO DEFAULT(NONE)&
!!$OMP PRIVATE(indt,dvd,dvu,uvd,uvu,f,g,x,flo,uuilon,duilon,ulouilopn)&
!!$OMP PRIVATE(cil,temp,wronk,i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmpl)&
!!$OMP PRIVATE(lmplm,lmx,lmxx,lp,lp1,lpl,m,mem,mems,mp,mu,n,nh,noded)&
!!$OMP PRIVATE(nodeu,nsym,na)&
!!$OMP SHARED(dimension,atoms,jspin,jsp,sphhar,enpara,td,ud,l_write,ci,vr,mpi,input)
       DO  n = 1,atoms%ntype
          na=sum(atoms%neq(:n-1))+1
          
          IF (l_write) WRITE (6,FMT=8000) n
          DO l = 0,atoms%lmax(n)
             CALL radfun(l,n,jspin,enpara%el0(l,n,jspin),vr(:,0,n),atoms,&
                  f(1,1,l),g(1,1,l),ud,nodeu,noded,wronk)
             IF (l_write) WRITE (6,FMT=8010) l,enpara%el0(l,n,jspin),ud%us(l,n,jspin),&
                  ud%dus(l,n,jspin),nodeu,ud%uds(l,n,jspin),ud%duds(l,n,jspin),noded,ud%ddn(l,n,jspin),wronk
             END DO
8000      FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',&
               /,t32,'radial function',t79,'energy derivative',/,t3,&
               'l',t8,'energy',t26,'value',t39,'derivative',t53,&
               'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,&
               'norm',t119,'wronskian')
8010      FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)
          !
          !--->   generate the extra wavefunctions for the local orbitals,
          !--->   if there are any.
          !
          IF (atoms%nlo(n).GE.1) THEN
             CALL radflo(atoms,n,jspin,enpara%ello0(1,1,jspin), vr(:,0,n), f,g,mpi,&
                  ud, uuilon,duilon,ulouilopn,flo)
          END IF
          
          nsym = atoms%ntypsy(na)
          nh = sphhar%nlh(nsym)
          !
          !--->    generate the irreducible integrals (u(l'):v(lamda,nu:u(l))
          !--->    for l' .ge. l, but only thos that will contribute
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
          
          td%tuu(:,n,jsp) = cmplx(0.0,0.0)
          td%tdd(:,n,jsp) = cmplx(0.0,0.0)
          td%tud(:,n,jsp) = cmplx(0.0,0.0)
          td%tdu(:,n,jsp) = cmplx(0.0,0.0)
          indt=0
          s=atoms%lmax(n)*(atoms%lmax(n)+2)+1
          !--->    generate the various t(l'm',lm) matrices for l'm'.ge.lm
          !--->    loop over l'm'
          DO lp = 0,atoms%lmax(n)
             lp1 = (lp* (lp+1))/2
             DO mp = -lp,lp
                lmp = lp* (lp+1) + mp
                lmpl = (lmp* (lmp+1))/2
                !--->    loop over lattice harmonics
                DO lh = 0, nh
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
          
          !Setup local hamiltonian
          DO lmp=0,atoms%lmax(n)*(atoms%lmax(n)+2)
             lp=FLOOR(SQRT(1.0*lmp))
             mp=lmp-lp*(lp+1)
             IF (lp>atoms%lmax(n).OR.ABS(mp)>lp) STOP "BUG"
             !--->             loop over l,m
             DO l = 0,atoms%lmax(n)
                DO m = -l,l
                   lm = l* (l+1) + m
                   in = td%ind(lmp,lm,n,jsp)
                   IF (in/=-9999) THEN
                      IF (in>=0) THEN
                         td%h_loc(lm,lmp,n,jsp) =CONJG(td%tuu(in,n,jsp))
                         td%h_loc(lm+s,lmp,n,jsp) =CONJG(td%tud(in,n,jsp))
                         td%h_loc(lm,lmp+s,n,jsp) =CONJG(td%tdu(in,n,jsp))
                         td%h_loc(lm+s,lmp+s,n,jsp) =CONJG(td%tdd(in,n,jsp))
                      ELSE
                         td%h_loc(lm,lmp,n,jsp) =td%tuu(-in,n,jsp)
                         td%h_loc(lm+s,lmp,n,jsp) =td%tdu(-in,n,jsp)
                         td%h_loc(lm,lmp+s,n,jsp) =td%tud(-in,n,jsp)
                         td%h_loc(lm+s,lmp+s,n,jsp) =td%tdd(-in,n,jsp)
                      END IF
                      !--->    update ax, bx
                      
                   END IF
                END DO
             END DO
          ENDDO
          !DO lp = 0,atoms%lmax(n)
          !   lp1 = (lp* (lp+1))/2+lp
          !   uvu(lp1,0)=uvu(lp1,0)+enpara%el0(lp,n,jsp)+td%e_shift
          !   uvd(lp1,0)=uvd(lp1,0)+0.5
          !   dvu(lp1,0)=dvu(lp1,0)+0.5
          !   dvd(lp1,0)=dvd(lp1,0)+(enpara%el0(lp,n,jsp)+td%e_shift)*ud%ddn(lp,n,jsp)
          !ENDDO
   
          !--->    include diagonal terms from muffin-tin hamiltonian
          DO lp = 0,atoms%lmax(n)
             DO mp = -lp,lp
                lmp = lp* (lp+1) + mp
                td%h_loc(lmp,lmp,n,jsp)=(td%e_shift+enpara%el0(lp,n,jsp))+td%h_loc(lmp,lmp,n,jsp)
                td%h_loc(lmp,lmp+s,n,jsp)=0.5+td%h_loc(lmp,lmp+s,n,jsp)
                td%h_loc(lmp+s,lmp,n,jsp)=0.5+td%h_loc(lmp+s,lmp,n,jsp)
                td%h_loc(lmp+s,lmp+s,n,jsp)=(enpara%el0(lp,n,jsp)+td%e_shift)*ud%ddn(lp,n,jsp)+td%h_loc(lmp+s,lmp+s,n,jsp)
             END DO
          END DO   

          !Perform cholesky decomposition
          CALL zpotrf("L",2*s,td%h_loc(:,:,n,jsp),size(td%h_loc,1),info)
          !Generate full matrix
          DO lm=0,size(td%h_loc,1)-1
             DO lmp=lm+1,size(td%h_loc,1)-1
                td%h_loc(lm,lmp,n,jsp)=0.0!td%h_loc(lm,lmp,n,jsp)
             ENDDO
          ENDDO
          IF (info.ne.0) THEN
             td%e_shift=td%e_shift*2.0
             print *,"Potential shift to small, increasing the value to:",td%e_shift
             if (td%e_shift>e_shift_max) call judft_error("Potential shift at maximum")
             OK=.false.
          ENDIF

          !
          !--->   set up the t-matrices for the local orbitals,
          !--->   if there are any
          IF (atoms%nlo(n).GE.1) THEN
             CALL tlo(atoms,sphhar,jspin,jsp,n,enpara,0,input,vr(1,0,n),&
                  na,flo,f,g,ud, uuilon,duilon,ulouilopn, td)
             
          ENDIF
          
       ENDDO
!!$OMP END PARALLEL DO
    enddo
  END SUBROUTINE tlmplm_cholesky

END MODULE m_tlmplm_cholesky
