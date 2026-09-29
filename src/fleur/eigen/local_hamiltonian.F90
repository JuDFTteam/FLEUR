!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_local_Hamiltonian
   USE m_judft
   IMPLICIT NONE
   PRIVATE
   PUBLIC:: local_ham, add_nonsph, extract_nonsph
  !*********************************************************************
  ! Sets up the k-independent muffin-tin Hamiltonian td%h in the unified
  ! radial basis (u, udot and LOs per l, see t_radfun and t_tlmplm):
  !   h = V_nonsph + H_sph + DFT+U/OPC (for LOs)
  ! and the non-spherical LAPW block td%h_loc_nonsph used by hsmt_nonsph,
  ! shifted to be positive definite and Cholesky decomposed.
  !*********************************************************************
CONTAINS
   SUBROUTINE local_ham(sphhar,atoms,sym,noco,nococonv,enpara,&
       fmpi,v,vx,inden,input,hub1inp,hub1data,td,alpha_hybrid,l_dfptmod,l_forces)
      !! l_dfptmod: no Cholesky decomposition
      !! l_forces:  LAPW part up to lmax and without DFT+U (forces add it separately)
      USE m_constants
      USE m_types

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

      INTEGER :: j1,j2,jsp,n,lh0
      COMPLEX :: one
      REAL, ALLOCATABLE :: vr(:,:)

      CALL timestart("local_hamiltonian")
      IF (input%secvar) CALL judft_error("Second variation is not supported",calledby="local_ham")
      CALL td%init(atoms,input%jspins,PRESENT(l_forces))

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
            !$OMP SHARED(atoms,sphhar,sym,enpara,j1,j2,jsp,v,vx,input,hub1inp,td,alpha_hybrid,one,noco)
            DO  n = 1,atoms%ntype
               IF (j1==j2.OR.noco%l_unrestrictMT(n).or.noco%l_constrained(n)) THEN
                  CALL nonsph_potential(atoms,enpara,v,vx,n,jsp,alpha_hybrid,vr,lh0)
                  CALL add_nonsph(td,n,atoms,sym,sphhar,input,hub1inp,vr,lh0,j1,j2,one)
               END IF
               CALL extract_nonsph(td,atoms,n,j1,j2)
               IF (jsp<3) CALL add_sph(td,n,jsp)
            ENDDO
            !$OMP end parallel do
            CALL add_ldaU(fmpi,inden,jsp,atoms,v,input,td,j1,j2,PRESENT(l_forces))
            ! For DFPT, do not decompose
            IF (jsp<3.AND..NOT.PRESENT(l_dfptmod)) CALL cholesky_decompose(td,atoms,jsp)
         END DO
      END DO

      IF (noco%l_soc.AND.noco%l_noco.AND..NOT.noco%l_ss) &
         CALL add_soc(fmpi,atoms,noco,nococonv,input,enpara,v,hub1inp,hub1data,td)

      CALL timestop("local_hamiltonian")
   END SUBROUTINE

   SUBROUTINE nonsph_potential(atoms,enpara,v,vx,n,iSpinV,alpha_hybrid,vr,lh0)
      !! potential entering the non-spherical integrals and the first lattice harmonic to use
      USE m_constants, ONLY: sfp_const
      USE m_types
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

   SUBROUTINE add_nonsph(td,n,atoms,sym,sphhar,input,hub1inp,vr,lh0,j1,j2,one)
      !! td%h(:,:,n,j1,j2) += one * sum_lh <r_i^{l'} Y_{l'm'}|V_lh|r_j^l Y_lm> for all slots i,j
      USE m_constants, ONLY: ImagUnit
      USE m_intgr, ONLY: intgr3
      USE m_gaunt, ONLY: gaunt1
      USE m_types
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_hub1inp),INTENT(IN)    :: hub1inp
      REAL,           INTENT(IN)    :: vr(:,0:)
      INTEGER,        INTENT(IN)    :: n,lh0,j1,j2
      COMPLEX,        INTENT(IN)    :: one

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

   LOGICAL FUNCTION l_nonsph_removed(atoms,input,hub1inp,n,l)
      USE m_types
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

   SUBROUTINE add_sph(td,n,jsp)
      !! spherical Hamiltonian, diagonal in lm
      USE m_types
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      INTEGER,        INTENT(IN)    :: n,jsp
      INTEGER :: l,m,nr
      REAL, ALLOCATABLE :: hs(:,:)
      DO l = 0,td%lrange(n)
         hs = td%radfun(n)%hsph(l,jsp)
         nr = td%radfun(n)%n_r(l)
         DO m = -l,l
            td%h(td%ind(:nr,l*(l+1)+m,n),td%ind(:nr,l*(l+1)+m,n),n,jsp,jsp) = &
               td%h(td%ind(:nr,l*(l+1)+m,n),td%ind(:nr,l*(l+1)+m,n),n,jsp,jsp) + hs
         END DO
      END DO
   END SUBROUTINE

   PURE INTEGER FUNCTION nonsph_size(atoms,n)
      USE m_types_atoms
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER,       INTENT(IN) :: n
      nonsph_size = atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+1
   END FUNCTION

   SUBROUTINE extract_nonsph(td,atoms,n,j1,j2)
      !! copy the u/udot block with l<=lnonsph of td%h into td%h_loc_nonsph
      USE m_types
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      INTEGER,        INTENT(IN)    :: n,j1,j2
      INTEGER, ALLOCATABLE :: idx(:)
      idx = nonsph_ind(td,atoms,n)
      td%h_loc_nonsph(:,:,n,j1,j2) = CMPLX(0.0,0.0)
      td%h_loc_nonsph(:SIZE(idx)-1,:SIZE(idx)-1,n,j1,j2) = td%h(idx,idx,n,j1,j2)
   END SUBROUTINE

   FUNCTION nonsph_ind(td,atoms,n) RESULT(idx)
      !! positions in td%h of the u and udot functions with l<=lnonsph
      USE m_types
      TYPE(t_tlmplm), INTENT(IN) :: td
      TYPE(t_atoms),  INTENT(IN) :: atoms
      INTEGER,        INTENT(IN) :: n
      INTEGER, ALLOCATABLE :: idx(:)
      idx = [td%ind(1,:nonsph_size(atoms,n)-1,n),td%ind(2,:nonsph_size(atoms,n)-1,n)]
   END FUNCTION

   FUNCTION nonsph_lm_ind(atoms,n,l) RESULT(idx)
      !! positions in h_loc_nonsph of u (1,:) and udot (2,:) for all m of l
      USE m_types_atoms
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER,       INTENT(IN) :: n,l
      INTEGER :: idx(2,-l:l),m
      idx(1,:) = [(l*(l+1)+m,m=-l,l)]
      idx(2,:) = idx(1,:)+nonsph_size(atoms,n)
   END FUNCTION

   SUBROUTINE add_lm_block(mat,idx,mm,r)
      !! mat(idx(i,m),idx(j,mp)) += mm(m,mp)*r(i,j)
      COMPLEX, INTENT(INOUT) :: mat(0:,0:)
      INTEGER, INTENT(IN)    :: idx(:,:)
      COMPLEX, INTENT(IN)    :: mm(:,:)
      REAL,    INTENT(IN)    :: r(:,:)
      INTEGER :: m,mp
      DO mp = 1,SIZE(mm,2)
         DO m = 1,SIZE(mm,1)
            mat(idx(:,m),idx(:,mp)) = mat(idx(:,m),idx(:,mp)) + mm(m,mp)*r
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE add_ldaU(fmpi,inden,jsp,atoms,v,input,td,j1,j2,l_forces)
      !! DFT+U, DFT+HIA and OPC; LOs get them only through h (u/udot parts)
      USE m_opc_setup
      USE m_types
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_potden),   INTENT(IN)    :: v,inden
      TYPE(t_tlmplm),   INTENT(INOUT) :: td
      INTEGER,          INTENT(IN)    :: jsp,j1,j2
      LOGICAL,          INTENT(IN)    :: l_forces

      INTEGER  :: i_u,i_opc,n,l,m
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
         CALL add_lm_block(td%h_loc_nonsph(:,:,n,j1,j2),nonsph_lm_ind(atoms,n,l),mm,td%radfun(n)%integral(1:2,1:2,l,j1,j2))
         IF (atoms%lda_u(i_u)%use_lo.AND..NOT.l_forces) &
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
         CALL add_lm_block(td%h_loc_nonsph(:,:,n,j1,j2),nonsph_lm_ind(atoms,n,l),mm,td%radfun(n)%integral(1:2,1:2,l,j1,j2))
         IF (.NOT.l_forces) &
            CALL add_lm_block(td%h(:,:,n,j1,j2),td%ind(1:2,l*l:l*l+2*l,n),mm,td%radfun(n)%integral(1:2,1:2,l,j1,j2))
      END DO
   END SUBROUTINE

   SUBROUTINE add_soc(fmpi,atoms,noco,nococonv,input,enpara,v,hub1inp,hub1data,td)
      ! Setup of the soc parameters for first-variation SOC and the resulting
      ! correction of the relativistic LOs' spherical Hamiltonian.
      USE m_constants
      USE m_types

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

   SUBROUTINE cholesky_decompose(td,atoms,jsp)
      !! shift the non-spherical LAPW block by e_shift*overlap until it is positive definite
      USE m_types
      TYPE(t_tlmplm), INTENT(INOUT) :: td
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      INTEGER,        INTENT(IN)    :: jsp

      REAL, PARAMETER :: e_shift_min=0.5
      REAL, PARAMETER :: e_shift_max=65.0

      INTEGER :: n,info,s,l,i,m
      COMPLEX,ALLOCATABLE :: mat(:,:),shift(:,:)

      td%e_shift(:,jsp) = e_shift_min
      DO n = 1,atoms%ntype
         s = nonsph_size(atoms,n)
         info = 1
         DO WHILE(info.NE.0)
            mat = td%h_loc_nonsph(:2*s-1,:2*s-1,n,jsp,jsp)
            DO l = 0,atoms%lnonsph(n)
               ALLOCATE(shift(2*l+1,2*l+1),source=CMPLX(0.0,0.0))
               DO m = 1,2*l+1
                  shift(m,m) = td%e_shift(n,jsp)
               END DO
               CALL add_lm_block(mat,nonsph_lm_ind(atoms,n,l),shift,td%radfun(n)%integral(1:2,1:2,l,jsp,jsp))
               DEALLOCATE(shift)
            END DO
            CALL zpotrf("L",2*s,mat,SIZE(mat,1),info)
            DO i = 1,2*s
               mat(:i-1,i) = 0.0
            END DO
            IF (info.NE.0) THEN
               td%e_shift(n,jsp) = td%e_shift(n,jsp)*2.0
               IF (td%e_shift(n,jsp)>e_shift_max) CALL judft_error("Potential shift at maximum")
            END IF
         END DO
         td%h_loc_nonsph(:2*s-1,:2*s-1,n,jsp,jsp) = mat
      END DO
   END SUBROUTINE

END MODULE m_local_Hamiltonian
