!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_local_Hamiltonian
   USE m_judft
   USE m_constants
   USE m_intgr, ONLY: intgr3
   USE m_gaunt, ONLY: gaunt1
   USE m_gradYlm, ONLY: Derivative
   USE m_types_atoms
   USE m_opc_setup
   USE m_types_enpara
   USE m_types_hub1data
   USE m_types_hub1inp
   USE m_types_input
   USE m_types_mpi
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_potden
   USE m_types_radfun
   USE m_types_sphhar
   USE m_types_sym
   USE m_types_tlmplm
   IMPLICIT NONE
   PRIVATE
   PUBLIC:: local_ham, add_nonsph, extract_nonsph
  !*********************************************************************
  ! Sets up the k-independent muffin-tin Hamiltonian td%h in the unified
  ! radial basis (u, udot and LOs per l, see t_radfun and t_tlmplm):
  !   h = V_nonsph + H_sph + DFT+U/OPC (for LOs)
  ! and the non-spherical LAPW block td%h_loc_nonsph used by hsmt_nonsph,
  ! shifted to be positive definite and Cholesky decomposed.
  ! For MetaGGAs h also contains -1/2 div(V_tau grad), and td%h_sph_extra
  ! the spherical terms for l>lnonsph that hsmt_sph adds.
  !*********************************************************************
CONTAINS
   SUBROUTINE local_ham(sphhar,atoms,sym,noco,nococonv,enpara,&
       fmpi,v,vx,inden,input,hub1inp,hub1data,td,alpha_hybrid,l_dfptmod,l_forces,vTau)
      !! l_dfptmod: no Cholesky decomposition
      !! l_forces:  LAPW part up to lmax and without DFT+U (forces add it separately)
      !! vTau:      MetaGGA potential dE_xc/dtau (collinear only)

      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_sphhar),   INTENT(IN)    :: sphhar
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_hub1inp),  INTENT(IN)    :: hub1inp
      TYPE(t_hub1data), INTENT(INOUT) :: hub1data
      TYPE(t_potden),   INTENT(IN)    :: v,vx,inden
      TYPE(t_tlmplm),   INTENT(INOUT) :: td
      REAL,             INTENT(IN)    :: alpha_hybrid
      LOGICAL, INTENT(IN),OPTIONAL    :: l_dfptmod, l_forces
      TYPE(t_potden),   INTENT(IN), OPTIONAL :: vTau

      INTEGER :: j1,j2,jsp,n,lh0
      COMPLEX :: one
      REAL, ALLOCATABLE :: vr(:,:)

      CALL timestart("local_hamiltonian")
      IF (input%secvar) CALL judft_error("Second variation is not supported",calledby="local_ham")
      CALL td%init(atoms,input%jspins,PRESENT(l_forces))
      td%l_sph_extra = PRESENT(vTau)

      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atoms,input,enpara,fmpi,v,hub1data,td)
      DO n = 1,atoms%ntype
         CALL td%radfun(n)%generate_radial_functions(atoms,input,enpara,fmpi,v,n,hub1data)
      END DO
      !$OMP END PARALLEL DO

      DO jsp=1,MERGE(4,input%jspins,any(noco%l_unrestrictMT).OR.any(noco%l_spinoffd_ldau).or.any(noco%l_constrained))
         DO j1=merge(jsp,1,jsp<3),merge(jsp,2,jsp<3)
            j2  = MERGE(j1,3-j1,jsp<3)
            one = MERGE(CMPLX(1.,0.),CMPLX(0.,1.),jsp<4)
            one = MERGE(CONJG(one),one,j1<j2)

            !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(vr,lh0)&
            !$OMP SHARED(atoms,sphhar,sym,enpara,j1,j2,jsp,v,vx,input,hub1inp,td,alpha_hybrid,one,noco,vTau)
            DO  n = 1,atoms%ntype
               IF (j1==j2.OR.noco%l_unrestrictMT(n).or.noco%l_constrained(n)) THEN
                  CALL nonsph_potential(atoms,enpara,v,vx,n,jsp,alpha_hybrid,vr,lh0)
                  IF (PRESENT(vTau).AND.jsp<3) THEN
                     ! lambda=0 carries V_tau and V_sph-enpara%vr: the radial functions solve the
                     ! auxiliary GGA potential (Eq. 8 of Doumont et al., PRB 105, 195138)
                     lh0 = 0
                     CALL add_nonsph(td,n,atoms,sym,sphhar,input,hub1inp,vr,lh0,j1,j2,one,vtau=vTau%mt(:,0:,n,jsp))
                     CALL add_sph_extra(td,n,atoms,sym,sphhar,vr(:,0),vTau%mt(:,0,n,jsp),jsp)
                  ELSE
                     CALL add_nonsph(td,n,atoms,sym,sphhar,input,hub1inp,vr,lh0,j1,j2,one)
                  END IF
               END IF
               CALL extract_nonsph(td,atoms,n,j1,j2)
               IF (jsp<3) CALL add_sph(td,n,jsp,input%l_useapw)
            ENDDO
            !$OMP end parallel do
            CALL add_ldaU(fmpi,inden,jsp,atoms,v,input,td,j1,j2,PRESENT(l_forces))
            ! For DFPT, do not decompose
            IF (jsp<3.AND..NOT.PRESENT(l_dfptmod)) CALL cholesky_decompose(fmpi,td,atoms,jsp)
         END DO
      END DO

      IF (noco%l_soc.AND.noco%l_noco.AND..NOT.noco%l_ss) &
         CALL add_soc(fmpi,atoms,noco,nococonv,input,enpara,v,hub1inp,hub1data,td)

      CALL timestop("local_hamiltonian")
   END SUBROUTINE

   SUBROUTINE nonsph_potential(atoms,enpara,v,vx,n,iSpinV,alpha_hybrid,vr,lh0)
      !! potential entering the non-spherical integrals and the first lattice harmonic to use
      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_potden), INTENT(IN) :: v,vx
      INTEGER,        INTENT(IN) :: n,iSpinV
      REAL,           INTENT(IN) :: alpha_hybrid
      REAL, ALLOCATABLE, INTENT(OUT) :: vr(:,:)
      INTEGER,        INTENT(OUT):: lh0

      INTEGER :: jri

      jri = atoms%jri(n)
      lh0 = MERGE(1,0,iSpinV<3.and.alpha_hybrid==0)
      IF (atoms%l_nonpolbas(n)) lh0 = 0
      ALLOCATE(vr(SIZE(v%mt,1),0:SIZE(v%mt,2)-1),source=0.0)
      vr(:jri,0:) = v%mt(:jri,0:,n,iSpinV)
      IF (iSpinV<3) THEN
         vr(:jri,0) = (v%mt(:jri,0,n,iSpinV)-enpara%vr(:jri,n,iSpinV))/atoms%rmsh(:jri,n)*sfp_const
         IF (alpha_hybrid.NE.0) vr(:jri,0:) = vr(:jri,0:)-alpha_hybrid*vx%mt(:jri,0:,n,iSpinV)
      END IF
   END SUBROUTINE

   SUBROUTINE add_nonsph(td,n,atoms,sym,sphhar,input,hub1inp,vr,lh0,j1,j2,one,vtau)
      !! td%h(:,:,n,j1,j2) += one * sum_lh <r_i^{l'} Y_{l'm'}|V_lh|r_j^l Y_lm> for all slots i,j
      !! vtau (MetaGGA, needs lh0=0): also the elements of -1/2 div(V_tau grad), see vtau_integrals
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_hub1inp),INTENT(IN)    :: hub1inp
      REAL,           INTENT(IN)    :: vr(:,0:)
      INTEGER,        INTENT(IN)    :: n,lh0,j1,j2
      COMPLEX,        INTENT(IN)    :: one
      REAL, OPTIONAL, INTENT(IN)    :: vtau(:,0:)

      REAL, ALLOCATABLE :: vint(:,:,:,:,:)
      COMPLEX :: cil
      INTEGER :: nsym,nh,lr,jri,lp,l,lh,lamda,i,j,mp,m,mem,mu,lmp,lm

      ASSOCIATE(rf => td%radfun(n))
      nsym = sym%ntypsy(atoms%firstAtom(n))
      nh   = sphhar%nlh(nsym)
      lr   = td%lrange(n)
      jri  = atoms%jri(n)

      ! radial integrals of all slot pairs allowed by the Gaunt selection rules
      ALLOCATE(vint(MAXVAL(rf%n_r),MAXVAL(rf%n_r),0:lr,0:lr,lh0:MAX(lh0,nh)),source=0.0)
      DO lh = lh0,nh
         lamda = sphhar%llh(lh,nsym)
         DO lp = 0,lr
            DO l = 0,lr
               IF (MOD(lamda+lp+l,2)==1 .OR. lamda<ABS(lp-l) .OR. lamda>lp+l) CYCLE
               DO j = 1,rf%n_r(l)
                  DO i = 1,rf%n_r(lp)
                     CALL intgr3((rf%r(:jri,1,i,lp,j1)*rf%r(:jri,1,j,l,j2)+rf%r(:jri,2,i,lp,j1)*rf%r(:jri,2,j,l,j2))*vr(:jri,lh),&
                                 atoms%rmsh(:,n),atoms%dx(n),jri,vint(i,j,lp,l,lh))
                  END DO
               END DO
            END DO
         END DO
      END DO
      ! DFT+U/HIA orbitals: non-spherical double counting removed from their u/udot block
      DO l = 0,lr
         IF (l_nonsph_removed(atoms,input,hub1inp,n,l)) vint(1:2,1:2,l,l,:) = 0.0
      END DO
      IF (PRESENT(vtau)) THEN
         IF (lh0/=0) CALL judft_bug("V_tau needs the lambda=0 lattice harmonic",calledby="add_nonsph")
         CALL vtau_integrals(atoms,sphhar,rf,n,nsym,lr,j1,j2,vtau,vint)
      END IF

      DO lp = 0,lr
         DO mp = -lp,lp
            lmp = lp*(lp+1)+mp
            DO lh = lh0,nh
               lamda = sphhar%llh(lh,nsym)
               DO mem = 1,sphhar%nmem(lh,nsym)
                  mu = sphhar%mlh(mem,lh,nsym)
                  m  = mp-mu
                  DO l = ABS(lp-lamda),MIN(lp+lamda,lr),2
                     IF (ABS(m)>l) CYCLE
                     lm  = l*(l+1)+m
                     cil = ImagUnit**(l-lp)*sphhar%clnu(mem,lh,nsym)*gaunt1(lp,lamda,l,mp,mu,m,atoms%lmaxd)
                     td%h(td%ind(:rf%n_r(lp),lmp,n),td%ind(:rf%n_r(l),lm,n),n,j1,j2) = &
                        td%h(td%ind(:rf%n_r(lp),lmp,n),td%ind(:rf%n_r(l),lm,n),n,j1,j2) + one*cil*vint(:rf%n_r(lp),:rf%n_r(l),lp,l,lh)
                  END DO
               END DO
            END DO
         END DO
      END DO
      END ASSOCIATE
   END SUBROUTINE

   SUBROUTINE vtau_integrals(atoms,sphhar,rf,n,nsym,lr,j1,j2,vtau,vint)
      !! vint(i,j,l',l,lh) += 1/2 int [D_i^{l'} D_j^l + a/r^2 r_i^{l'} r_j^l] V_tau,lh dr, the radial part of
      !! 1/2 <grad(r_i^{l'} Y_{l'm'})|V_tau,lh Y_lh|grad(r_j^l Y_lm)>. D = dr/dr - r/r for the large and
      !! small component; a = [l(l+1)+l'(l'+1)-lambda(lambda+1)]/2 comes from the angular gradients.
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_radfun), INTENT(IN)    :: rf
      INTEGER,        INTENT(IN)    :: n,nsym,lr,j1,j2
      REAL,           INTENT(IN)    :: vtau(:,0:)
      REAL,           INTENT(INOUT) :: vint(:,:,0:,0:,0:)

      REAL, ALLOCATABLE :: d(:,:,:,:,:)
      REAL    :: a,x
      INTEGER :: jri,lh,lamda,lp,l,i,j,c,s,js(2)

      jri = atoms%jri(n)
      js  = [j1,j2]
      ALLOCATE(d(jri,2,MAXVAL(rf%n_r),0:lr,2))
      DO s = 1,2
         DO l = 0,lr
            DO i = 1,rf%n_r(l)
               DO c = 1,2
                  CALL Derivative(rf%r(:jri,c,i,l,js(s)),n,atoms,d(:,c,i,l,s))
                  d(:,c,i,l,s) = d(:,c,i,l,s)-rf%r(:jri,c,i,l,js(s))/atoms%rmsh(:jri,n)
               END DO
            END DO
         END DO
      END DO

      DO lh = 0,sphhar%nlh(nsym)
         lamda = sphhar%llh(lh,nsym)
         DO lp = 0,lr
            DO l = 0,lr
               IF (MOD(lamda+lp+l,2)==1 .OR. lamda<ABS(lp-l) .OR. lamda>lp+l) CYCLE
               a = 0.5*REAL(l*(l+1)+lp*(lp+1)-lamda*(lamda+1))
               DO j = 1,rf%n_r(l)
                  DO i = 1,rf%n_r(lp)
                     CALL intgr3((d(:,1,i,lp,1)*d(:,1,j,l,2)+d(:,2,i,lp,1)*d(:,2,j,l,2) &
                                 +a*(rf%r(:jri,1,i,lp,j1)*rf%r(:jri,1,j,l,j2)+rf%r(:jri,2,i,lp,j1)*rf%r(:jri,2,j,l,j2)) &
                                 /atoms%rmsh(:jri,n)**2)*vtau(:jri,lh),atoms%rmsh(:,n),atoms%dx(n),jri,x)
                     vint(i,j,lp,l,lh) = vint(i,j,lp,l,lh)+0.5*x
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE add_sph_extra(td,n,atoms,sym,sphhar,vr0,vtau0,jsp)
      !! MetaGGA, l>td%lrange(n): lambda=0 elements of V_sph-enpara%vr (vr0, as in nonsph_potential) and of
      !! -1/2 div(V_tau grad) for u and udot, which hsmt_sph adds to its energy-parameter based setup.
      !! Order in td%h_sph_extra: <u|u>, <u|udot>, <udot|u>, <udot|udot>.
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      REAL,           INTENT(IN)    :: vr0(:),vtau0(:)
      INTEGER,        INTENT(IN)    :: n,jsp

      REAL    :: d(atoms%jri(n),2,2),e(2,2),c0
      INTEGER :: jri,nsym,l,i,j,c

      jri  = atoms%jri(n)
      nsym = sym%ntypsy(atoms%firstAtom(n))
      ASSOCIATE(rf => td%radfun(n))
      DO l = td%lrange(n)+1,atoms%lmax(n)
         c0 = REAL(sphhar%clnu(1,0,nsym))*gaunt1(l,0,l,0,0,0,atoms%lmaxd)
         DO i = 1,2
            DO c = 1,2
               CALL Derivative(rf%r(:jri,c,i,l,jsp),n,atoms,d(:,c,i))
               d(:,c,i) = d(:,c,i)-rf%r(:jri,c,i,l,jsp)/atoms%rmsh(:jri,n)
            END DO
         END DO
         DO j = 1,2
            DO i = 1,2
               CALL intgr3((rf%r(:jri,1,i,l,jsp)*rf%r(:jri,1,j,l,jsp)+rf%r(:jri,2,i,l,jsp)*rf%r(:jri,2,j,l,jsp)) &
                           *(vr0(:jri)+0.5*REAL(l*(l+1))*vtau0(:jri)/atoms%rmsh(:jri,n)**2) &
                           +0.5*(d(:,1,i)*d(:,1,j)+d(:,2,i)*d(:,2,j))*vtau0(:jri),atoms%rmsh(:,n),atoms%dx(n),jri,e(i,j))
            END DO
         END DO
         td%h_sph_extra(:,l,n,jsp) = c0*[e(1,1),e(1,2),e(2,1),e(2,2)]
      END DO
      END ASSOCIATE
   END SUBROUTINE

   LOGICAL FUNCTION l_nonsph_removed(atoms,input,hub1inp,n,l)
      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_input),  INTENT(IN) :: input
      TYPE(t_hub1inp),INTENT(IN) :: hub1inp
      INTEGER,        INTENT(IN) :: n,l
      INTEGER :: i
      l_nonsph_removed = .FALSE.
      DO i = 1,atoms%n_u+atoms%n_hia
         IF (atoms%lda_u(i)%atomType/=n .OR. atoms%lda_u(i)%l/=l) CYCLE
         IF (i<=atoms%n_u) l_nonsph_removed = l_nonsph_removed .OR. input%ldauNonsphDC
         IF (i> atoms%n_u) l_nonsph_removed = l_nonsph_removed .OR. hub1inp%l_nonsphDC
      END DO
   END FUNCTION

   SUBROUTINE add_sph(td,n,jsp,l_useapw)
      !! spherical Hamiltonian, diagonal in lm; with APW including the kinetic surface
      !! term, which hsmt_sph adds itself for the LAPW part (h_loc_nonsph is extracted before)
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      INTEGER,        INTENT(IN)    :: n,jsp
      LOGICAL,        INTENT(IN)    :: l_useapw
      INTEGER :: l,m,nr
      REAL, ALLOCATABLE :: hs(:,:)
      DO l = 0,td%lrange(n)
         hs = td%radfun(n)%hsph(l,jsp)
         IF (l_useapw) hs = hs + td%radfun(n)%hsurf(l,jsp)
         nr = td%radfun(n)%n_r(l)
         DO m = -l,l
            td%h(td%ind(:nr,l*(l+1)+m,n),td%ind(:nr,l*(l+1)+m,n),n,jsp,jsp) = &
               td%h(td%ind(:nr,l*(l+1)+m,n),td%ind(:nr,l*(l+1)+m,n),n,jsp,jsp) + hs
         END DO
      END DO
   END SUBROUTINE

   PURE INTEGER FUNCTION nonsph_size(atoms,n)
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER,       INTENT(IN) :: n
      nonsph_size = atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+1
   END FUNCTION

   SUBROUTINE extract_nonsph(td,atoms,n,j1,j2)
      !! copy the block of u and LAPW udot with l<=lnonsph of td%h into td%h_loc_nonsph
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      INTEGER,        INTENT(IN)    :: n,j1,j2
      INTEGER, ALLOCATABLE :: idx(:)
      idx = nonsph_ind(td,atoms,n)
      td%h_loc_nonsph(:,:,n,j1,j2) = CMPLX(0.0,0.0)
      td%h_loc_nonsph(:SIZE(idx)-1,:SIZE(idx)-1,n,j1,j2) = td%h(idx,idx,n,j1,j2)
   END SUBROUTINE

   FUNCTION nonsph_ind(td,atoms,n) RESULT(idx)
      !! positions in td%h of the u and the LAPW udot functions with l<=lnonsph,
      !! in the order of the matching coefficients of hsmt_ab
      TYPE(t_tlmplm), INTENT(IN) :: td
      TYPE(t_atoms),  INTENT(IN) :: atoms
      INTEGER,        INTENT(IN) :: n
      INTEGER, ALLOCATABLE :: idx(:)
      INTEGER :: l
      idx = td%ind(1,:nonsph_size(atoms,n)-1,n)
      DO l = 0,atoms%lnonsph(n)
         IF (.NOT.atoms%l_apw(l,n)) idx = [idx,td%ind(2,l*l:l*l+2*l,n)]
      END DO
   END FUNCTION

   FUNCTION nonsph_lm_ind(atoms,n,l) RESULT(idx)
      !! positions in h_loc_nonsph of u (1,:) and, for LAPW channels, udot (2,:) for all m of l
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER,       INTENT(IN) :: n,l
      INTEGER, ALLOCATABLE :: idx(:,:)
      INTEGER :: boff(0:atoms%lnonsph(n))
      INTEGER :: m
      boff(0:atoms%lnonsph(n)) = atoms%udot_rows(atoms%lnonsph(n),n)
      ALLOCATE(idx(MERGE(1,2,boff(l)<0),2*l+1))
      idx(1,:) = [(l*(l+1)+m,m=-l,l)]
      IF (boff(l)>=0) idx(2,:) = [(boff(l)+l+m,m=-l,l)]
   END FUNCTION

   SUBROUTINE add_nonsph_lm_block(mat,atoms,n,l,mm,r)
      !! h_loc_nonsph-shaped mat += mm(m,mp)*r(i,j) for the u (and LAPW udot) functions of l
      COMPLEX,       INTENT(INOUT) :: mat(0:,0:)
      TYPE(t_atoms), INTENT(IN)    :: atoms
      INTEGER,       INTENT(IN)    :: n,l
      COMPLEX,       INTENT(IN)    :: mm(:,:)
      REAL,          INTENT(IN)    :: r(:,:)
      INTEGER, ALLOCATABLE :: idx(:,:)
      idx = nonsph_lm_ind(atoms,n,l)
      CALL add_lm_block(mat,idx,mm,r(:SIZE(idx,1),:SIZE(idx,1)))
   END SUBROUTINE add_nonsph_lm_block

   SUBROUTINE add_lm_block(mat,idx,mm,r)
      !! mat(idx(i,m),idx(j,mp)) += mm(m,mp)*r(i,j)
      COMPLEX, INTENT(INOUT) :: mat(0:,0:)
      INTEGER, INTENT(IN)    :: idx(:,:)
      COMPLEX, INTENT(IN)    :: mm(:,:)
      REAL,    INTENT(IN)    :: r(:,:)
      INTEGER :: m,mp,i,j
      DO mp = 1,SIZE(mm,2)
         DO m = 1,SIZE(mm,1)
            DO j = 1,SIZE(idx,1)
               DO i = 1,SIZE(idx,1)
                  mat(idx(i,m),idx(j,mp)) = mat(idx(i,m),idx(j,mp)) + mm(m,mp)*r(i,j)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE add_ldaU(fmpi,inden,jsp,atoms,v,input,td,j1,j2,l_forces)
      !! DFT+U, DFT+HIA and OPC; LOs get them only through h (u/udot parts)
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_potden),   INTENT(IN)    :: v,inden
      TYPE(t_tlmplm),   INTENT(INOUT) :: td
      INTEGER,          INTENT(IN)    :: jsp,j1,j2
      LOGICAL,          INTENT(IN)    :: l_forces

      INTEGER  :: i_u,i_opc,n,l,m
      LOGICAL  :: l_apwlo
      COMPLEX, ALLOCATABLE :: mm(:,:)
      REAL, ALLOCATABLE :: opc_corrections(:)

      DO i_u = 1,atoms%n_u+atoms%n_hia
         n = atoms%lda_u(i_u)%atomType
         l = atoms%lda_u(i_u)%l
         IF (j1==j2) THEN
            mm = v%mmpMat(-l:l,-l:l,i_u,jsp)
         ELSE IF (j1>j2) THEN
            mm = v%mmpMat(-l:l,-l:l,i_u,3)
         ELSE
            mm = CONJG(TRANSPOSE(v%mmpMat(-l:l,-l:l,i_u,3)))
         END IF
         CALL add_nonsph_lm_block(td%h_loc_nonsph(:,:,n,j1,j2),atoms,n,l,mm,td%radfun(n)%integral(1:2,1:2,l,j1,j2))
         ! An APW LO carries most of the weight of its channel and gets U through td%h;
         ! spin off-diagonal blocks only in their last pass (jsp=4), after extract_nonsph
         l_apwlo = ANY(atoms%l_dulo(:atoms%nlo(n),n).AND.atoms%llo(:atoms%nlo(n),n)==l)
         IF ((atoms%lda_u(i_u)%use_lo.OR.(l_apwlo.AND.(j1==j2.OR.jsp==4))).AND..NOT.l_forces) &
            CALL add_lm_block(td%h(:,:,n,j1,j2),td%ind(1:2,l*l:l*l+2*l,n),mm,td%radfun(n)%integral(1:2,1:2,l,j1,j2))
      END DO

      IF (atoms%n_opc==0 .OR. j1/=j2) RETURN
      CALL opc_setup(input,atoms,fmpi,v,inden,jsp,opc_corrections)
      DO i_opc = 1,atoms%n_opc
         n = atoms%lda_opc(i_opc)%atomType
         l = atoms%lda_opc(i_opc)%l
         IF (ALLOCATED(mm)) DEALLOCATE(mm)
         ALLOCATE(mm(2*l+1,2*l+1),source=CMPLX(0.0,0.0))
         DO m = -l,l
            mm(m+l+1,m+l+1) = opc_corrections(i_opc)*m
         END DO
         CALL add_nonsph_lm_block(td%h_loc_nonsph(:,:,n,j1,j2),atoms,n,l,mm,td%radfun(n)%integral(1:2,1:2,l,j1,j2))
         IF (.NOT.l_forces) &
            CALL add_lm_block(td%h(:,:,n,j1,j2),td%ind(1:2,l*l:l*l+2*l,n),mm,td%radfun(n)%integral(1:2,1:2,l,j1,j2))
      END DO
   END SUBROUTINE

   SUBROUTINE add_soc(fmpi,atoms,noco,nococonv,input,enpara,v,hub1inp,hub1data,td)
      ! Setup of the soc parameters for first-variation SOC and the resulting
      ! correction of the relativistic LOs' spherical Hamiltonian.

      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_potden),   INTENT(IN)    :: v
      TYPE(t_hub1inp),  INTENT(IN)    :: hub1inp
      TYPE(t_hub1data), INTENT(INOUT) :: hub1data
      TYPE(t_tlmplm),   INTENT(INOUT) :: td

      INTEGER :: n,l,m,jsp,lo,i_hia,i

      ! Radial SOC matrix rsoc%rso used by hsmt_soc_offdiag
      IF (.NOT.ALLOCATED(td%rsoc%rso)) CALL td%rsoc%init(atoms)
      CALL td%rsoc%rad_matrix(atoms,noco,nococonv,input,fmpi,enpara,v)
      ! Hubbard-1 SOC parameter xi from <u|V_SO|u>
      IF (fmpi%irank==0) THEN
         DO i_hia = 1, atoms%n_hia
            IF (hub1inp%l_soc_given(i_hia)) CYCLE
            n = atoms%lda_u(atoms%n_u+i_hia)%atomType
            l = atoms%lda_u(atoms%n_u+i_hia)%l
            hub1data%xi(i_hia) = 2.0*td%rsoc%rso(1,1,n,l,1,1)*hartree_to_ev_const
         END DO
      END IF

      ! relLO: its own spherical element is the Dirac eigenvalue, which already contains
      ! the SOC of the j=l-1/2 branch. H_sph must carry the SOC-free value, as the
      ! first-variational SOC is added by hsmt_soc:
      !     <relLO|H_sph|relLO> = epsilon + (l+1)*I_so,   I_so = rsoc%rso(slot,slot,n,l)
      DO n = 1, atoms%ntype
         DO lo = 1, atoms%nlo(n)
            IF (.NOT.atoms%l_relLO(lo,n)) CYCLE
            l = atoms%llo(lo,n)
            DO jsp = 1, input%jspins
               DO m = -l, l
                  i = td%ind(atoms%slot_of_lo(lo,n),l*(l+1)+m,n)
                  td%h(i,i,n,jsp,jsp) = td%h(i,i,n,jsp,jsp) + REAL(l+1)*td%rsoc%rso(atoms%slot_of_lo(lo,n),atoms%slot_of_lo(lo,n),n,l,jsp,jsp)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE cholesky_decompose(fmpi,td,atoms,jsp)
      !! shift the non-spherical LAPW block by e_shift*overlap until it is positive definite
      TYPE(t_mpi),    INTENT(IN)    :: fmpi
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      INTEGER,        INTENT(IN)    :: jsp

      REAL, PARAMETER :: e_shift_min=0.5
      REAL, PARAMETER :: e_shift_max=65.0

      INTEGER :: n,info,s,l,i,m
      COMPLEX,ALLOCATABLE :: mat(:,:),shift(:,:)

      td%e_shift(:,jsp) = e_shift_min
      DO n = 1,atoms%ntype
         s = atoms%num_ab_rows(atoms%lnonsph(n),n)
         info = 1
         DO WHILE(info.NE.0)
            mat = td%h_loc_nonsph(:s-1,:s-1,n,jsp,jsp)
            DO l = 0,atoms%lnonsph(n)
               ALLOCATE(shift(2*l+1,2*l+1),source=CMPLX(0.0,0.0))
               DO m = 1,2*l+1
                  shift(m,m) = td%e_shift(n,jsp)
               END DO
               CALL add_nonsph_lm_block(mat,atoms,n,l,shift,td%radfun(n)%integral(1:2,1:2,l,jsp,jsp))
               DEALLOCATE(shift)
            END DO
            CALL zpotrf("L",s,mat,SIZE(mat,1),info)
            DO i = 1,s
               mat(:i-1,i) = 0.0
            END DO
            IF (info.NE.0) THEN
               td%e_shift(n,jsp) = td%e_shift(n,jsp)*2.0
               IF (td%e_shift(n,jsp)>e_shift_max) THEN
                  CALL cholesky_failure_report(fmpi,td,atoms,n,jsp,info)
                  CALL judft_error("Potential shift at maximum",calledby="cholesky_decompose",&
                                   hint="Local Hamiltonian not positive definite, see diagnostics above")
               END IF
            END IF
         END DO
         td%h_loc_nonsph(:s-1,:s-1,n,jsp,jsp) = mat
      END DO
   END SUBROUTINE

   SUBROUTINE cholesky_failure_report(fmpi,td,atoms,n,jsp,info)
      !! diagnostics of the unshifted non-spherical LAPW block written before aborting
      TYPE(t_mpi),    INTENT(IN) :: fmpi
      TYPE(t_tlmplm), INTENT(IN) :: td
      TYPE(t_atoms),  INTENT(IN) :: atoms
      INTEGER,        INTENT(IN) :: n,jsp,info

      INTEGER :: s,ns,k,l,m,ierr,nbad
      INTEGER :: boff(0:atoms%lnonsph(n))
      REAL, ALLOCATABLE    :: eig(:),rwork(:)
      COMPLEX, ALLOCATABLE :: h(:,:),work(:)

      s = nonsph_size(atoms,n)
      ns = atoms%num_ab_rows(atoms%lnonsph(n),n)
      boff(0:atoms%lnonsph(n)) = atoms%udot_rows(atoms%lnonsph(n),n)
      h = td%h_loc_nonsph(:ns-1,:ns-1,n,jsp,jsp)
      WRITE(*,'(a,i0,a,i0,a,i0,a,i0)') "Rank ",fmpi%irank,": Cholesky decomposition of local Hamiltonian failed for atom type ",&
         n,", spin ",jsp,", lnonsph ",atoms%lnonsph(n)
      WRITE(*,'(a,f8.3,a,i0)') "  last shift (Htr): ",td%e_shift(n,jsp),", zpotrf info: ",info
      IF (info>0) THEN
         IF (info<=s) THEN
            k = info-1
            l = INT(SQRT(REAL(k)+0.5))
         ELSE
            l = MAXLOC(boff,1,MASK=boff<=info-1)-1
            k = info-1-boff(l)+l*l
         END IF
         m = k-l*(l+1)
         WRITE(*,'(3a,i0,a,i0,a,2es14.5)') "  failing basis function: ",MERGE("u   ","udot",info<=s)," l=",l," m=",m,&
            ", unshifted diagonal element: ",h(info,info)
      END IF
      nbad = COUNT(.NOT.(ABS(h)<=HUGE(1.0)))
      WRITE(*,'(a,i0)') "  NaN/Inf matrix elements: ",nbad
      IF (nbad==0) THEN
         WRITE(*,'(a,es14.5)') "  max |H-H^H|: ",MAXVAL(ABS(h-CONJG(TRANSPOSE(h))))
         ALLOCATE(eig(ns),rwork(3*ns),work(2*ns))
         CALL zheev("N","L",ns,h,ns,eig,work,SIZE(work),rwork,ierr)
         IF (ierr==0) WRITE(*,'(a,4es14.5)') "  lowest eigenvalues (unshifted): ",eig(:MIN(4,ns))
      END IF
      WRITE(*,'(a)') "  radial overlaps  l   <u|u>          <u|udot>       <udot|udot>"
      DO l = 0,atoms%lnonsph(n)
         WRITE(*,'(15x,i3,3es15.5)') l,td%radfun(n)%integral(1,1,l,jsp,jsp),td%radfun(n)%integral(1,2,l,jsp,jsp),&
            td%radfun(n)%integral(2,2,l,jsp,jsp)
      END DO
   END SUBROUTINE

END MODULE m_local_Hamiltonian
