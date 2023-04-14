MODULE m_tlmplm_cholesky
   USE m_judft

   IMPLICIT NONE
  !*********************************************************************
  !     sets up the local Hamiltonian, i.e. the Hamiltonian in the
  !     l',m',l,m,u- basis which is independent from k!
  !     shifts this local Hamiltonian to make it positive definite
  !     and does a cholesky decomposition
  !*********************************************************************
CONTAINS
   SUBROUTINE tlmplm_cholesky(sphhar,atoms,sym,noco,nococonv,enpara,&
       jspin,fmpi,v,vx,inden,input,hub1inp,hub1data,td,ud,alpha_hybrid,l_dfptmod)
      !! Should probably be called tlmplm_postprocess or something, as there is a
      !! lot more happening than only the Cholesky decomposition.

      ! Add onto the t_{L'L}^{\mu} matrices from tlmplm some contributions from
      ! DFT+U, DFT+HIA, DFT+OPC, constraint fields and the diagonal E_{l} terms
      ! etc. In the diagonal case, Cholesky decompose the nonspherical part of
      ! the Hamiltonian by shifting the diagonal part upwards until the matrix
      ! is positive-definite.

      USE m_tlmplm
      USE m_types
      USE m_radovlp
      USE m_opc_setup

      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_sphhar),   INTENT(IN)    :: sphhar
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_hub1inp),  INTENT(IN)    :: hub1inp
      TYPE(t_hub1data), INTENT(IN)    :: hub1data
      TYPE(t_potden),   INTENT(IN)    :: v,vx,inden
      TYPE(t_tlmplm),   INTENT(INOUT) :: td
      TYPE(t_usdus),    INTENT(INOUT) :: ud

      ! Scalar Arguments
      INTEGER, INTENT(IN) :: jspin  ! Spin index for v
      REAL,    INTENT(IN) :: alpha_hybrid

      LOGICAL, INTENT(IN),OPTIONAL :: l_dfptmod

      ! Local Scalars
      INTEGER :: i,l,lm,lmp,lp,info,in,jsp,j1,j2
      INTEGER :: mp,n,m,s,i_u,jmin,jmax,i_opc
      LOGICAL :: OK, isRoot, l_call_tlmplm, l_dfpt
      COMPLEX :: one

      ! Local Arrays
      REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)
      REAL, ALLOCATABLE :: opc_corrections(:)

      REAL, PARAMETER :: e_shift_min=0.5
      REAL, PARAMETER :: e_shift_max=65.0

      l_dfpt = PRESENT(l_dfptmod)

      jsp = jspin
      IF (jsp<3) THEN
         ! e_shift for Cholesky decomposition
         jmin=jsp;jmax=jsp
         td%e_shift(:,jsp)=e_shift_min
      ELSE
         jmin=1;jmax=2
         IF(atoms%n_u+atoms%n_hia>0) THEN
            !Calculate overlap integrals for the off-diagonal LDA+U contributions
            ALLOCATE(uun21(0:atoms%lmaxd,atoms%ntype),source=0.0)
            ALLOCATE(udn21(0:atoms%lmaxd,atoms%ntype),source=0.0)
            ALLOCATE(dun21(0:atoms%lmaxd,atoms%ntype),source=0.0)
            ALLOCATE(ddn21(0:atoms%lmaxd,atoms%ntype),source=0.0)
            CALL rad_ovlp(atoms,ud,input,hub1data,v%mt,enpara%el0, uun21,udn21,dun21,ddn21)
         END IF
      END IF

      isRoot = fmpi%is_root()

      IF (atoms%n_opc > 0 .and. jsp<3) THEN
         CALL opc_setup(input,atoms,fmpi,v,inden,jsp,opc_corrections)
      END IF

      DO j1=jmin,jmax
         j2  = MERGE(j1,3-j1,jsp<3)
         one = MERGE(CMPLX(1.,0.),CMPLX(0.,1.),jsp<4)
         one = MERGE(CONJG(one),one,j1<j2)

         l_call_tlmplm = j1.EQ.j2.OR.any(noco%l_unrestrictMT)

         !$OMP PARALLEL DO DEFAULT(NONE)&
         !$OMP PRIVATE(i,l,lm,lmp)&
         !$OMP PRIVATE(lp,m,mp,n)&
         !$OMP PRIVATE(OK,s,in,info,i_u,i_opc)&
         !$OMP SHARED(one,nococonv,atoms,jspin,jsp,sym,sphhar,enpara,td,ud,v,vx,alpha_hybrid,isRoot,l_call_tlmplm)&
         !$OMP SHARED(fmpi,input,hub1inp,hub1data,uun21,udn21,dun21,ddn21,opc_corrections,j1,j2,l_dfpt)
         DO  n = 1,atoms%ntype
            IF (l_call_tlmplm) THEN

               CALL tlmplm(n,sphhar,atoms,sym,enpara,nococonv,j1,j2,jsp,fmpi,v,vx,input,hub1inp,hub1data,td,ud,alpha_hybrid,one,l_dfpt)
            END IF
            OK = .FALSE.
            cholesky_loop:DO WHILE(.NOT.OK)
               OK = .TRUE.
               s = td%h_loc2_nonsph(n)
               ! Set up local hamiltonian
               td%h_loc_nonsph(0:s-1,0:s-1,n,j1,j2)    = td%h_loc(0:s-1,0:s-1,n,j1,j2)
               td%h_loc_nonsph(s:s+s-1,0:s-1,n,j1,j2)  = td%h_loc(td%h_loc2(n):s+td%h_loc2(n)-1,0:s-1,n,j1,j2)
               td%h_loc_nonsph(0:s-1,s:s+s-1,n,j1,j2)  = td%h_loc(0:s-1,td%h_loc2(n):s+td%h_loc2(n)-1,n,j1,j2)
               td%h_loc_nonsph(s:s+s-1,s:s+s-1,n,j1,j2)= td%h_loc(td%h_loc2(n):s+td%h_loc2(n)-1,td%h_loc2(n):s+td%h_loc2(n)-1,n,j1,j2)

               ! Include contribution from LDA+U and LDA+HIA (latter are behind LDA+U contributions)
               DO i_u=1,atoms%n_u+atoms%n_hia
                  IF (n.NE.atoms%lda_u(i_u)%atomType) CYCLE
                  ! Found a "U" for this atom type
                  l  = atoms%lda_u(i_u)%l
                  lp = atoms%lda_u(i_u)%l
                  DO m = -l,l
                     lm = l* (l+1) + m
                     DO mp = -lp,lp
                        lmp = lp*(lp+1) + mp
                        IF (j1==j2) THEN
                           td%h_loc_nonsph(lm  ,lmp  ,n,j1,j2) = td%h_loc_nonsph(lm  ,lmp  ,n,j1,j2) + v%mmpMat(m,mp,i_u,jsp)
                           td%h_loc_nonsph(lm+s,lmp+s,n,j1,j2) = td%h_loc_nonsph(lm+s,lmp+s,n,j1,j2) + v%mmpMat(m,mp,i_u,jsp) * ud%ddn(lp,n,jsp)
                        ELSE IF(j1>j2) THEN
                           td%h_loc_nonsph(lm  ,lmp  ,n,j1,j2) = td%h_loc_nonsph(lm  ,lmp  ,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * uun21(l,n)
                           td%h_loc_nonsph(lm  ,lmp+s,n,j1,j2) = td%h_loc_nonsph(lm  ,lmp+s,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * udn21(l,n)
                           td%h_loc_nonsph(lm+s,lmp  ,n,j1,j2) = td%h_loc_nonsph(lm+s,lmp  ,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * dun21(l,n)
                           td%h_loc_nonsph(lm+s,lmp+s,n,j1,j2) = td%h_loc_nonsph(lm+s,lmp+s,n,j1,j2) + v%mmpMat(m,mp,i_u,3) * ddn21(l,n)
                        ELSE
                           ! For this part of the Hamiltonian we need to perform Hermitian conjugation on mmpMat
                           td%h_loc_nonsph(lm  ,lmp  ,n,j1,j2) = td%h_loc_nonsph(lm  ,lmp  ,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * uun21(l,n)
                           td%h_loc_nonsph(lm  ,lmp+s,n,j1,j2) = td%h_loc_nonsph(lm  ,lmp+s,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * udn21(l,n)
                           td%h_loc_nonsph(lm+s,lmp  ,n,j1,j2) = td%h_loc_nonsph(lm+s,lmp  ,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * dun21(l,n)
                           td%h_loc_nonsph(lm+s,lmp+s,n,j1,j2) = td%h_loc_nonsph(lm+s,lmp+s,n,j1,j2) + conjg(v%mmpMat(mp,m,i_u,3)) * ddn21(l,n)
                        END IF
                     END DO
                  END DO
               END DO

               DO i_opc=1,atoms%n_opc
                  IF (n.NE.atoms%lda_opc(i_opc)%atomType) CYCLE
                  ! Found an "OPC" for this atom type
                  l=atoms%lda_opc(i_opc)%l
                  DO m = -l,l
                     lm = l*(l+1) + m
                     td%h_loc_nonsph(lm  ,lm  ,n,j1,j2) = td%h_loc_nonsph(lm  ,lm  ,n,j1,j2) + opc_corrections(i_opc) * m
                     td%h_loc_nonsph(lm+s,lm+s,n,j1,j2) = td%h_loc_nonsph(lm+s,lm+s,n,j1,j2) + opc_corrections(i_opc) * m * ud%ddn(l,n,jsp)
                  END DO
               END DO

               ! Create Cholesky decomposition of local hamiltonian
               ! For DFPT, do not decompose!
               IF (jsp<3.AND..NOT.l_dfpt) THEN
                  ! Add shift onto the diagonal terms to make matrix positive definite
                  DO lp = 0,atoms%lnonsph(n)
                     DO mp = -lp,lp
                        lmp = lp* (lp+1) + mp
                        td%h_loc_nonsph(lmp,lmp,n,j1,j2)=td%e_shift(n,jsp)+td%h_loc_nonsph(lmp,lmp,n,j1,j2)
                        td%h_loc_nonsph(lmp+s,lmp+s,n,j1,j2)=td%e_shift(n,jsp)*ud%ddn(lp,n,jsp)+td%h_loc_nonsph(lmp+s,lmp+s,n,j1,j2)
                     END DO
                  END DO
                  IF (lmp+1.NE.s) CALL judft_error("BUG in tlmplm_cholesky")

                  ! Perform cholesky decomposition
                  info=0
                  CALL zpotrf("L",2*s,td%h_loc_nonsph(:,:,n,j1,j2),SIZE(td%h_loc_nonsph,1),info)

                  ! Set upper part to zero
                  DO l=0,2*s-1
                     DO lp=0,l-1
                        td%h_loc_nonsph(lp,l,n,j1,j2)=0.0
                     END DO
                  END DO

                  IF (info.NE.0) THEN
                     td%e_shift(n,jsp)=td%e_shift(n,jsp)*2.0
                     IF (isRoot) THEN
                        PRINT *,"Potential shift for atom type ",n," is too small. Increasing the value to:",td%e_shift(n,jsp)
                     END IF
                     IF (td%e_shift(n,jsp)>e_shift_max) THEN
                        CALL judft_error("Potential shift at maximum")
                     END IF
                     OK = .FALSE.
                  END IF
               END IF
            END DO cholesky_loop

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
                     ! For DFPT we need a non-symmetrized local Hamiltonian
                     !IF (l_dfpt) THEN
                     !   td%h_loc(lm,lm+s,n,jsp,jsp) = td%h_loc(lm,lm+s,n,jsp,jsp) + 0.5 !+1.0
                     !   td%h_loc(lm+s,lm,n,jsp,jsp) = td%h_loc(lm+s,lm,n,jsp,jsp) - 0.5 !+0.0
                     !END IF
                  END DO
               END DO
            END IF
         END DO
         !$OMP END PARALLEL DO
         IF (any(noco%l_constrained)) CALL tlmplm_constrained(atoms,v,enpara,input,hub1data,ud,nococonv,td)
      END DO

   END SUBROUTINE tlmplm_cholesky

   SUBROUTINE tlmplm_constrained(atoms,v,enpara,input,hub1data,ud,nococonv,td)
      USE m_radovlp
      USE m_types

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

END MODULE m_tlmplm_cholesky
