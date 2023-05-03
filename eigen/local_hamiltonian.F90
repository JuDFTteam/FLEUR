MODULE m_local_Hamiltonian
   USE m_judft
   IMPLICIT NONE
   PRIVATE
   PUBLIC:: local_ham
  !*********************************************************************
  !     sets up the local Hamiltonian, i.e. the Hamiltonian in the
  !     l',m',l,m,u- basis which is independent from k!
  !     shifts this local Hamiltonian to make it positive definite
  !     and does a cholesky decomposition
  !*********************************************************************
CONTAINS
   SUBROUTINE local_ham(sphhar,atoms,sym,noco,nococonv,enpara,&
       fmpi,v,vx,inden,input,hub1inp,hub1data,td,ud,alpha_hybrid,l_dfptmod)
      !! Should probably be called tlmplm_postprocess or something, as there is a
      !! lot more happening than only the Cholesky decomposition.

      ! Add onto the t_{L'L}^{\mu} matrices from tlmplm some contributions from
      ! DFT+U, DFT+HIA, DFT+OPC, constraint fields and the diagonal E_{l} terms
      ! etc. In the diagonal case, Cholesky decompose the nonspherical part of
      ! the Hamiltonian by shifting the diagonal part upwards until the matrix
      ! is positive-definite.
      USE m_spnorb
      USE m_tlmplm
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
      TYPE(t_usdus),    INTENT(INOUT) :: ud

      ! Scalar Arguments
      
      REAL,    INTENT(IN) :: alpha_hybrid

      LOGICAL, INTENT(IN),OPTIONAL :: l_dfptmod

      ! Local Scalars
      INTEGER :: l,lm,j1,j2,jsp
      INTEGER :: n,m,s
      COMPLEX :: one

      CALL timestart("local_hamiltonian")
      CALL td%init(atoms,input%jspins,(noco%l_noco.AND.noco%l_soc.AND..NOT.noco%l_ss).OR.any(noco%l_constrained))

      DO jsp=1,MERGE(4,input%jspins,any(noco%l_unrestrictMT).OR.any(noco%l_spinoffd_ldau))

         DO j1=merge(jsp,1,jsp<3),merge(jsp,2,jsp<3)
            j2  = MERGE(j1,3-j1,jsp<3)
            one = MERGE(CMPLX(1.,0.),CMPLX(0.,1.),jsp<4)
            one = MERGE(CONJG(one),one,j1<j2)

            !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(l,m,lm,s)&
            !$OMP SHARED(atoms,sphhar,sym,enpara,nococonv,j1,j2,jsp,fmpi,v,vx,input,hub1inp,hub1data,td,ud,alpha_hybrid,one,l_dfptmod,noco)
            DO  n = 1,atoms%ntype
               IF (j1==j2.OR.noco%l_unrestrictMT(n)) THEN
                  CALL tlmplm(n,sphhar,atoms,sym,enpara,nococonv,j1,j2,jsp,fmpi,v,vx,input,hub1inp,hub1data,td,ud,alpha_hybrid,one,PRESENT(l_dfptmod))
               END IF
               !Copy local hamiltonian for non_spherical setup
               call restrict_to_lnonsph(td%h_loc(:,:,n,j1,j2),td%h_loc2(n),td%h_loc2_nonsph(n),td%h_loc_nonsph(:,:,n,j1,j2))
               ! Now add diagonal contributions to the matrices:
               IF (jsp<3) THEN
                  DO l = 0,atoms%lmax(n)
                     DO  m = -l,l
                        lm = l*(l+1) + m
                        s = td%h_loc2(n)
                        td%h_loc(lm,lm,n,jsp,jsp)     = td%h_loc(lm,lm,n,jsp,jsp)     + enpara%el0(l,n,jsp)
                        td%h_loc(lm,lm+s,n,jsp,jsp)   = td%h_loc(lm,lm+s,n,jsp,jsp)   + 0.5 ! Symmetrized from 1.0
                        td%h_loc(lm+s,lm,n,jsp,jsp)   = td%h_loc(lm+s,lm,n,jsp,jsp)   + 0.5 ! Symmetrized from 0.0
                        td%h_loc(lm+s,lm+s,n,jsp,jsp) = td%h_loc(lm+s,lm+s,n,jsp,jsp) + enpara%el0(l,n,jsp)*ud%ddn(l,n,jsp)
                     END DO
                  END DO
               END IF
               ! Store Matrices needed for LOs 
               if (atoms%nlotot>0) call restrict_to_lnonsph(td%h_loc(:,:,n,j1,j2),td%h_loc2(n),td%h_loc2_nonsph(n),td%h_loc_LO(:,:,n,j1,j2))
            ENDDO
            !$OMP end parallel do
            !Add LDA+U
            call add_ldaU(fmpi,inden,jsp,atoms,v,ud,input,hub1data,enpara,td%h_loc_nonsph,td%h_loc2_nonsph,j1,j2,.false.)
            !Add LDA+U also to LO part   
            if (atoms%nlotot>0) call add_ldaU(fmpi,inden,jsp,atoms,v,ud,input,hub1data,enpara,td%h_loc_LO,td%h_loc2_nonsph,j1,j2,.true.)
            ! Create Cholesky decomposition of local hamiltonian
            ! For DFPT, do not decompose!
            IF (jsp<3.AND..NOT.PRESENT(l_dfptmod)) THEN
               call cholesky_decompose(td%h_loc_nonsph(:,:,:,j1,j2),td%e_shift(:,jsp),atoms,ud,jsp)
            endif

            IF (any(noco%l_constrained)) CALL tlmplm_constrained(atoms,v,enpara,input,hub1data,ud,nococonv,td)
         END DO
      END DO

      !Setup of soc parameters for first-variation SOC
      IF (noco%l_soc.AND.noco%l_noco.AND..NOT.noco%l_ss) THEN
         CALL spnorb(atoms,noco,nococonv,input,fmpi,enpara,v%mt,ud,td%rsoc,.FALSE.,hub1inp,hub1data)
      END IF
      CALL timestop("local_hamiltonian")
      

   END SUBROUTINE 

   SUBROUTINE tlmplm_constrained(atoms,v,enpara,input,hub1data,ud,nococonv,td)
      USE m_radovlp
      USE m_types
      USE m_tlmplm
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_potden),   INTENT(IN)    :: v
      TYPE(t_tlmplm),   INTENT(INOUT) :: td
      TYPE(t_usdus),    INTENT(INOUT) :: ud
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_hub1data), INTENT(IN)    :: hub1data

      REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)
 
      COMPLEX :: c
      INTEGER :: n,l,s

      ALLOCATE(uun21(0:atoms%lmaxd,atoms%ntype),udn21(0:atoms%lmaxd,atoms%ntype), &
             & dun21(0:atoms%lmaxd,atoms%ntype),ddn21(0:atoms%lmaxd,atoms%ntype))

      CALL rad_ovlp(atoms,ud,input,hub1data,v%mt,enpara%el0,uun21,udn21,dun21,ddn21)

      DO  n = 1,atoms%ntype
         ! In a constraint calculation, we have to calculate the local spin
         ! off-diagonal contributions

         s=atoms%lnonsph(n)+1
         ! First ispin=1,jspin=2 case
         DO l=0, atoms%lnonsph(n)
            c = (-0.5)*CMPLX(nococonv%b_con(1,n),nococonv%b_con(2,n))
            td%h_off(l  ,l  ,n,1,2)     =td%h_off(l  ,l  ,n,1,2) + uun21(l,n)*c
            td%h_off(l  ,l+s,n,1,2)     =td%h_off(l  ,l+s,n,1,2) + udn21(l,n)*c
            td%h_off(l+s,l  ,n,1,2)     =td%h_off(l+s,l  ,n,1,2) + dun21(l,n)*c
            td%h_off(l+s,l+s,n,1,2)     =td%h_off(l+s,l+s,n,1,2) + ddn21(l,n)*c
         END DO

         ! Then ispin=2,jspin=1 case
         DO l = 0, atoms%lnonsph(n)
            c = (-0.5)*CMPLX(nococonv%b_con(1,n),-nococonv%b_con(2,n))
            td%h_off(l  ,l  ,n,2,1)     =td%h_off(l  ,l  ,n,2,1) + uun21(l,n)*c
            td%h_off(l  ,l+s,n,2,1)     =td%h_off(l  ,l+s,n,2,1) + udn21(l,n)*c
            td%h_off(l+s,l  ,n,2,1)     =td%h_off(l+s,l  ,n,2,1) + dun21(l,n)*c
            td%h_off(l+s,l+s,n,2,1)     =td%h_off(l+s,l+s,n,2,1) + ddn21(l,n)*c
         END DO
      END DO
   END SUBROUTINE tlmplm_constrained

   SUBROUTINE add_ldaU(fmpi,inden,jsp,atoms,v,ud,input,hub1data,enpara,mat,mat_half,j1,j2,l_lomatrix)
               ! Include contribution from LDA+U and LDA+HIA (latter are behind LDA+U contributions)
      USE m_radovlp
      USE m_opc_setup
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_hub1data), INTENT(IN)    :: hub1data
      TYPE(t_usdus),    INTENT(INOUT) :: ud
      TYPE(t_potden),   INTENT(IN)    :: v,inden
      INTEGER,INTENT(IN)              :: jsp,j1,j2
      COMPLEX,INTENT(INOUT)           :: mat(:,:,:,:,:)
      INTEGER,INTENT(IN)              :: mat_half(:)
      LOGICAL,INTENT(IN)              :: l_lomatrix

      INTEGER  :: i_u,n,s,l,lp,m,lm,mp,lmp,i_opc
      REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)
      REAL, ALLOCATABLE :: opc_corrections(:)

      IF(atoms%n_u+atoms%n_hia+atoms%n_opc==0) return !No LDA+U

      IF (j1.ne.j2) THEN
         !Calculate overlap integrals for the off-diagonal LDA+U contributions
         ALLOCATE(uun21(0:atoms%lmaxd,atoms%ntype),source=0.0)
         ALLOCATE(udn21(0:atoms%lmaxd,atoms%ntype),source=0.0)
         ALLOCATE(dun21(0:atoms%lmaxd,atoms%ntype),source=0.0)
         ALLOCATE(ddn21(0:atoms%lmaxd,atoms%ntype),source=0.0)
         CALL rad_ovlp(atoms,ud,input,hub1data,v%mt,enpara%el0, uun21,udn21,dun21,ddn21)
      ELSE 
         IF (atoms%n_opc > 0) THEN
            CALL opc_setup(input,atoms,fmpi,v,inden,jsp,opc_corrections)
         END IF
      
      END IF
   
      !!$OMP PARALLEL DO DEFAULT(NONE) private(n,s,l,lp,m,lm,mp,lmp)&
      !!$OMP shared(atoms,mat_half,l_lomatrix,v,mat,uun21,udn21,dun21,ddn21,j1,j2,ud,jsp)
      DO i_u=1,atoms%n_u+atoms%n_hia
         if (l_lomatrix.and..not.atoms%lda_u(i_u)%use_lo) cycle
         n=atoms%lda_u(i_u)%atomType
         s=mat_half(n)
         ! Found a "U" for this atom type
         l  = atoms%lda_u(i_u)%l
         lp = atoms%lda_u(i_u)%l
         
         DO m = -l,l
            lm = l* (l+1) + m +1 !indexing from 1
            DO mp = -lp,lp
               lmp = lp*(lp+1) + mp +1 !indexing from 1
               IF (j1==j2) THEN
                  mat(lm  ,lmp  ,n,j1,j2) = mat(lm  ,lmp  ,n,j1,j2) + v%mmpMat(m,mp,i_u,jsp)
                  mat(lm+s,lmp+s,n,j1,j2) = mat(lm+s,lmp+s,n,j1,j2) + v%mmpMat(m,mp,i_u,jsp) * ud%ddn(lp,n,jsp)
               ELSE IF(j1>j2) THEN
                  mat(lm  ,lmp  ,n,j1,j2) = mat(lm  ,lmp  ,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * uun21(l,n)
                  mat(lm  ,lmp+s,n,j1,j2) = mat(lm  ,lmp+s,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * udn21(l,n)
                  mat(lm+s,lmp  ,n,j1,j2) = mat(lm+s,lmp  ,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * dun21(l,n)
                  mat(lm+s,lmp+s,n,j1,j2) = mat(lm+s,lmp+s,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * ddn21(l,n)
               ELSE
                  ! For this part of the Hamiltonian we need to perform Hermitian conjugation on mmpMat
                  mat(lm  ,lmp  ,n,j1,j2) = mat(lm  ,lmp  ,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * uun21(l,n)
                  mat(lm  ,lmp+s,n,j1,j2) = mat(lm  ,lmp+s,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * udn21(l,n)
                  mat(lm+s,lmp  ,n,j1,j2) = mat(lm+s,lmp  ,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * dun21(l,n)
                  mat(lm+s,lmp+s,n,j1,j2) = mat(lm+s,lmp+s,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * ddn21(l,n)
               END IF
            END DO
         END DO
      END DO
      !!$OMP end parallel do
      DO i_opc=1,atoms%n_opc
         n=atoms%lda_opc(i_opc)%atomType
         s=mat_half(n)
         ! Found an "OPC" for this atom type
         l=atoms%lda_opc(i_opc)%l
         DO m = -l,l
            lm = l*(l+1) + m +1 !indexing from 1
            mat(lm  ,lm  ,n,j1,j2) = mat(lm  ,lm  ,n,j1,j2) + opc_corrections(i_opc) * m
            mat(lm+s,lm+s,n,j1,j2) = mat(lm+s,lm+s,n,j1,j2) + opc_corrections(i_opc) * m * ud%ddn(l,n,jsp)
         END DO
      END DO
    END SUBROUTINE

    subroutine cholesky_decompose(matrix,e_shift,atoms,ud,jsp)
    USE m_types
    TYPE(t_atoms),    INTENT(IN)    :: atoms
    TYPE(t_usdus),INTENT(IN)        :: ud
    COMPLEX,INTENT(INOUT)           :: matrix(0:,0:,:)
    REAL, INTENT(OUT)               :: e_shift(:)
    INTEGER,INTENT(IN)              :: jsp

    REAL, PARAMETER :: e_shift_min=0.5
    REAL, PARAMETER :: e_shift_max=65.0

    INTEGER :: n,info,s,l,lp,lmp,mp
    COMPLEX,ALLOCATABLE :: mat(:,:)

    e_shift=e_shift_min

    DO n=1,atoms%ntype
      s=atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+1
      info=1
      cholesky_loop:DO WHILE(info.ne.0)
         mat=matrix(0:,0:,n)
         !Mat is now using a lower bound of 1!!
         ! Add shift onto the diagonal terms to make matrix positive definite
         DO lp = 0,atoms%lnonsph(n)
               DO mp = -lp,lp
               lmp = lp* (lp+1) + mp +1
               mat(lmp,lmp)=e_shift(n)+mat(lmp,lmp)
               mat(lmp+s,lmp+s)=e_shift(n)*ud%ddn(lp,n,jsp)+mat(lmp+s,lmp+s)
               END DO
         END DO
         IF (lmp.NE.s) CALL judft_error("BUG in local_Hamiltonian:cholesky")

         ! Perform cholesky decomposition
         CALL zpotrf("L",2*s,mat(:,:),SIZE(mat,1),info)

         ! Set upper part to zero
         DO l=1,2*s
               DO lp=1,l-1
               mat(lp,l)=0.0
               END DO
         END DO

         IF (info.NE.0) THEN
               e_shift(n)=e_shift(n)*2.0
               IF (e_shift(n)>e_shift_max) THEN
               CALL judft_error("Potential shift at maximum")
               END IF
         END IF
      END DO cholesky_loop
      matrix(0:,0:,n)=mat
   ENDDO   
   END SUBROUTINE


    subroutine restrict_to_lnonsph(mat,s2,s,mat_nonsph)
        COMPLEX,INTENT(IN)   :: mat(0:,0:)
        INTEGER,INTENT(IN)   :: s,s2
        COMPLEX,INTENT(OUT)  :: mat_nonsph(0:,0:)
        ! Set up local hamiltonian
        mat_nonsph(0:s-1,0:s-1)    = mat(0:s-1,0:s-1)
        mat_nonsph(s:s+s-1,0:s-1)  = mat(s2:s+s2-1,0:s-1)
        mat_nonsph(0:s-1,s:s+s-1)  = mat(0:s-1,s2:s+s2-1)
        mat_nonsph(s:s+s-1,s:s+s-1)= mat(s2:s+s2-1,s2:s+s2-1)
    end subroutine
        
END MODULE m_local_Hamiltonian
