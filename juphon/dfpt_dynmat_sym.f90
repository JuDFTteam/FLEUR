MODULE m_dfpt_dynmat_sym
   USE m_juDFT
   USE m_types
   USE m_constants

   IMPLICIT NONE
CONTAINS
   SUBROUTINE ft_dyn(atoms, qpts, sym, amat, dyn_mat_q, dyn_mat_r, dyn_mat_q_full)
      !! Transforms the dynamical matrices for a set of q vectors in the
      !! irreducible Brillouin zone onto the full set of q vector in the BZ
      !! and subsequently transforms it to real space (lattice vector grid), to
      !! calculate the mass-normalized Force Constant Matrix.
      type(t_atoms), INTENT(IN)   :: atoms
      type(t_kpts),  INTENT(IN)   :: qpts
      type(t_sym),   INTENT(IN)   :: sym
      REAL,          INTENT(IN)   :: amat(3,3)
      COMPLEX,       INTENT(IN)  :: dyn_mat_q(:,:,:) ! (nqpt,dyn_dim,dyn_dim)
      COMPLEX,       INTENT(OUT) :: dyn_mat_r(:,:,:) ! (nqptf,dyn_dim,dyn_dim)

      COMPLEX, ALLOCATABLE, INTENT(OUT) :: dyn_mat_q_full(:,:,:)

      INTEGER :: iq, dyn_dim, iqfull
      INTEGER :: isym, r_lim(2,3)
      INTEGER :: j1, j2, j3
      REAL    :: q_full(3), q_full_BZ(3)

      INTEGER, ALLOCATABLE :: qvec_to_index(:,:,:)
      COMPLEX, ALLOCATABLE :: dyn_mat_qsym(:,:)

      dyn_dim = 3*atoms%nat

      allocate(qvec_to_index(0:qpts%nkpt3(1)-1,0:qpts%nkpt3(2)-1,0:qpts%nkpt3(3)-1))
      allocate(dyn_mat_qsym(dyn_dim,dyn_dim))
      allocate(dyn_mat_q_full(qpts%nkptf,dyn_dim,dyn_dim))

      qvec_to_index = 0

      ! Create an array that maps the q coordinates to the index of their q vector
      ! in the full BZ
      DO j1 = 0, qpts%nkpt3(1)-1
         q_full(1) = j1*1.0/qpts%nkpt3(1)
         DO j2 = 0, qpts%nkpt3(2)-1
            q_full(2) = j2*1.0/qpts%nkpt3(2)
            DO j3 = 0, qpts%nkpt3(3)-1
               q_full(3) = j3*1.0/qpts%nkpt3(3)
               DO iq = 1, qpts%nkptf
                  IF (norm2(q_full-qpts%bkf(:,iq))<1e-8) qvec_to_index(j1, j2, j3) = iq
               END DO
            END DO
         END DO
      END DO

      r_lim(2,:) = qpts%nkpt3(:)/2
      r_lim(1,:) = r_lim(2,:) - qpts%nkpt3(:) + 1

      dyn_mat_r(:,:,:) = cmplx(0.0,0.0)

      DO j1 = 0, qpts%nkpt3(1)-1
         q_full(1) = j1*1.0/qpts%nkpt3(1)
         DO j2 = 0, qpts%nkpt3(2)-1
            q_full(2) = j2*1.0/qpts%nkpt3(2)
            DO j3 = 0, qpts%nkpt3(3)-1
               q_full(3) = j3*1.0/qpts%nkpt3(3)
               ! Get q vector index and that of its representative in the irreducible wedge
               iqfull = qvec_to_index(j1, j2, j3)
               iq = qpts%bkp(iqfull)
               ! Fold vector back to 1st BZ if necessary
               q_full_BZ(:)=q_full(:)
               q_full_BZ = qpts%to_first_bz(q_full)

               isym = sym%invtab(qpts%bksym(iqfull))
               dyn_mat_qsym(:,:) = cmplx(0.0,0.0)
               CALL rotate_dynmat(atoms,sym,isym,amat,qpts%bk(:,iq),dyn_mat_q(iq,:,:),dyn_mat_qsym)
               dyn_mat_qsym(:,:) = dyn_mat_qsym(:,:)
               dyn_mat_q_full(iqfull,:,:) = dyn_mat_qsym

               ! Perform the actual FT onto the lattice vector grid
               CALL ft_dyn_direct(r_lim,1,q_full,dyn_mat_qsym,dyn_mat_r)
            END DO
         END DO
      END DO
      dyn_mat_r(:,:,:)=dyn_mat_r(:,:,:)/qpts%nkptf
   END SUBROUTINE

   SUBROUTINE ft_dyn_direct(ft_lim,isn,bqpt,dyn_mat_q,dyn_mat_r)
      INTEGER, INTENT(IN) :: ft_lim(2,3), isn
      REAL,    INTENT(IN) :: bqpt(3)

      COMPLEX :: dyn_mat_q(:,:)
      COMPLEX :: dyn_mat_r(:,:,:)

      INTEGER :: iGrid, ix, iy, iz, iout
      REAL    :: phas
      COMPLEX :: phase_fac
      iGrid=0
      DO iz=ft_lim(1,3),ft_lim(2,3)
         DO iy=ft_lim(1,2),ft_lim(2,2)
            DO ix=ft_lim(1,1),ft_lim(2,1)
               iGrid = iGrid+1
               phas=isn*tpi_const*(bqpt(1)*ix+bqpt(2)*iy+bqpt(3)*iz)
               phase_fac=cmplx(cos(phas),sin(phas))
               IF (isn==1) THEN
                  dyn_mat_r(iGrid,:,:) = dyn_mat_r(iGrid,:,:) + phase_fac*dyn_mat_q(:,:)
               ELSE IF (isn==-1) THEN
                  dyn_mat_q(:,:)    = dyn_mat_q(:,:)    + phase_fac*dyn_mat_r(iGrid,:,:)
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE rotate_dynmat(atoms,sym,isym,amat,bqpt,dyn,dyn_mat_qsym)
      !! Applies a symmetry operation to the dynamical matrix of an IBZ q vector
      !! to find the matrix of its mapped q vector in the full BZ. This is done
      !! by using the symmetry relation of the FCM when a symmetry operation
      !! \(\underline{B}\) maps atoms \(\beta',\alpha'\) onto \(\beta,\alpha\)
      !! in Cartesian coordinates 
      !! (\(\underline{B}\boldsymbol{\tau}_{\beta'}=\boldsymbol{\tau}_{\beta}\))
      !! $$\underline{\Phi}_{\alpha'+\boldsymbol{R},\beta'}=\underline{B}\underline{\Phi}_{\alpha+\underline{B}\boldsymbol{R},\beta}\underline{B}^{-1}$$
      !! Resulting (with the definition of the DM as the mass-scaled Fourier Transform
      !! of the FCM) in the corresponding relation:
      !! $$\underline{D}_{\alpha',\beta'}(\boldsymbol{q})=p(\alpha,\beta)*\underline{B}\underline{D}_{\alpha,\beta}(\boldsymbol{q}_{\mathrm{rep}})\underline{B}^{-1}$$
      !! with a phase factor
      !! $$f(\alpha,\beta)=exp(ix),x=\boldsymbol{q}\cdot(\boldsymbol{\tau}_{\beta'}-\boldsymbol{\tau}_{\alpha'})-\boldsymbol{q}_{\mathrm{red}}\cdot(\boldsymbol{\tau}_{\beta}-\boldsymbol{\tau}_{\alpha})$$
      !! which can be written as a product
      !! $$f(\alpha,\beta) = f(\alpha)f^{*}(\beta), f(\alpha)=exp(i(\boldsymbol{q}_{\mathrm{red}}\cdot\boldsymbol{\tau}_{\alpha}-\boldsymbol{q}\cdot\boldsymbol{\tau}_{\alpha'})).$$
      !! The real space rotation is related to the rotation matrix of the symmetry operation by the Bravais matrix of the system
      !! $$\underline{B}=\underline{A}\underline{S}^{-1}\underline{A}^{-1}$$
      USE m_inv3
      type(t_atoms), INTENT(IN)    :: atoms
      type(t_sym),   INTENT(IN)    :: sym
      INTEGER,       INTENT(IN)    :: isym
      REAL,          INTENT(IN)    :: amat(3,3)
      REAL,          INTENT(IN)    :: bqpt(3)
      COMPLEX,       INTENT(IN)    :: dyn(:,:)
      COMPLEX,       INTENT(INOUT) :: dyn_mat_qsym(:,:)

      INTEGER :: iAtom, jAtom
      INTEGER :: iAlpha, iBeta, iAlpha_map, iBeta_map
      REAL    :: mrot(3,3), invmrot(3,3), invamat(3,3), q_full(3), phas, det
      COMPLEX :: brot(3,3), temp_mat_1(3,3), temp_mat_2(3,3)
      COMPLEX :: phase_fac

      INTEGER :: map(atoms%nat)
      COMPLEX :: phase_map(atoms%nat)

      mrot = sym%mrot(:,:,isym)
      invmrot = sym%mrot(:,:,sym%invtab(isym))

      CALL inv3(amat,invamat,det)
      temp_mat_1 = MATMUL(invmrot,invamat)
      brot = MATMUL(amat,temp_mat_1)
      ! Find the q vector in the full BZ
      q_full = MATMUL(mrot,bqpt)

      ! Calculate the array of phases
      DO iAtom = 1, atoms%nat
         jAtom = sym%mapped_atom(isym,iAtom)
         map(iAtom)=jAtom
         ! Calculate the phase factor f(alpha)
         phas=tpi_const*(dot_product(bqpt(:),atoms%taual(:,iAtom))-dot_product(q_full(:),atoms%taual(:,jAtom)))
         phase_fac=cmplx(cos(phas),sin(phas))
         phase_map(iAtom)=phase_fac
      END DO
      ! Transform the dynamical matrix from the representative atom and q vector to the unfolded ones
      DO iAtom=1, atoms%nat
         iAlpha = 3*(iAtom-1)
         iAlpha_map = 3*(map(iAtom)-1)
         DO jAtom=1, atoms%nat
            iBeta = 3*(jAtom-1)
            iBeta_map = 3*(map(jAtom)-1)
            temp_mat_1 = dyn(iAlpha+1:iAlpha+3,iBeta+1:iBeta+3)
            temp_mat_2 = MATMUL(brot,temp_mat_1)
            temp_mat_1 = MATMUL(temp_mat_2,TRANSPOSE(brot))
            phase_fac=phase_map(iAtom)*conjg(phase_map(jAtom))

            dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) &
          = dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) + phase_fac*temp_mat_1
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE ift_dyn(atoms,qpts,sym,amat,bqpt,dyn_mat_r,dyn_mat_q)
      !! Transforms the dynamical matrix on a real space lattice vector grid
      !! (--> mass-normalized FCM) back onto a specific q vector provided as
      !! input (bqpt) by the inverse Fourier Transformation as compared to
      !! SUBROUTINE ft_dyn.
      type(t_atoms), INTENT(IN) :: atoms
      type(t_kpts),  INTENT(IN) :: qpts
      type(t_sym),   INTENT(IN) :: sym
      REAL,          INTENT(IN) :: amat(3,3)
      REAL,          INTENT(IN) :: bqpt(3)
      COMPLEX,       INTENT(IN) :: dyn_mat_r(:,:,:)

      COMPLEX, ALLOCATABLE, INTENT(OUT) :: dyn_mat_q(:,:)

      INTEGER :: dyn_dim, q_lim(2,3)

      dyn_dim = 3*atoms%nat

      allocate(dyn_mat_q(dyn_dim,dyn_dim))
      dyn_mat_q(:,:) = cmplx(0.0,0.0)

      q_lim(2,:) = qpts%nkpt3(:)/2
      q_lim(1,:) = q_lim(2,:) - qpts%nkpt3(:) + 1

      CALL ft_dyn_direct(q_lim,-1,bqpt,dyn_mat_q,dyn_mat_r)
   END SUBROUTINE

END MODULE m_dfpt_dynmat_sym