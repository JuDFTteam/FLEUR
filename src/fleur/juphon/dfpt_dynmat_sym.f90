MODULE m_dfpt_dynmat_sym
   USE m_juDFT
   USE m_types
   USE m_constants
   !USE m_npy

   IMPLICIT NONE
CONTAINS
   SUBROUTINE ft_dyn(atoms, qpts, sym, ft_lim, amat, l_WSinterpol, dyn_mat_q, dyn_mat_r, fft_grid, dyn_mat_q_full)
      !! Transforms the dynamical matrices for a set of q vectors in the
      !! irreducible Brillouin zone onto the full set of q vector in the BZ
      !! and subsequently transforms it to real space (lattice vector grid), to
      !! calculate the mass-normalized Force Constant Matrix.
      type(t_atoms), intent(in)       :: atoms
      type(t_kpts),  intent(in)       :: qpts
      type(t_sym),   intent(in)       :: sym
      integer,       intent(in)       :: ft_lim(2,3)
      real,          intent(in)       :: amat(3,3)
      logical, intent(in)             :: l_WSinterpol 
      complex,       intent(in)       :: dyn_mat_q(:,:,:) ! (nqpt,dyn_dim,dyn_dim)
      complex,intent(out)             :: dyn_mat_r(0:,0:,0:,:,:) ! (n1,n2,n3,dyn_dim,dyn_dim)
      complex, intent(out), allocatable :: fft_grid(:,:,:) ! (nqptf,dyn_dim,dyn_dim)
      complex, allocatable, intent(out) :: dyn_mat_q_full(:,:,:)

      integer :: mrot(3,3),invmrot(3,3)
      logical :: l_inv

      integer :: iq, dyn_dim, iqfull
      integer :: isym
      integer :: iz, iy, ix, iGrid, nx, ny, nz 
      real    :: q_full(3), trans(3)
      complex, allocatable :: dyn_mat_qsym(:,:)

      dyn_dim = 3*atoms%nat

      allocate(dyn_mat_qsym(dyn_dim,dyn_dim))
      allocate(dyn_mat_q_full(qpts%nkptf,dyn_dim,dyn_dim))
      allocate(fft_grid(qpts%nkptf,dyn_dim,dyn_dim))
      fft_grid(:,:,:) = cmplx(0.0,0.0)

      do iqfull = 1, qpts%nkptf
         ! Get q vector index and that of its representative in the irreducible wedge
         iq = qpts%bkp(iqfull)
         ! Fold vector back to 1st BZ if necessary
         q_full = qpts%bkf(:,iqfull)
         isym = qpts%bksym(iqfull)
         call sym%get_sym_operation_int_coord(isym,mrot,invmrot,trans,l_inv)
         if (.not. all(trans == 0 )) call juDFT_error("dynMat interpolation with non symmorphic symmetries is currently not supported. & 
                                                      Please redo the calculation of the IBZ and interpolation with a symmorphic group.",calledby="dfpt_dynmat_sym.f90")
         if (l_inv) isym = isym - sym%nop ! the corresponding symmetry operation 
         dyn_mat_qsym(:,:) = cmplx(0.0,0.0)
         call rotate_dynmat(atoms,sym,isym,mrot,invmrot,l_inv,amat,qpts%bk(:,iq),dyn_mat_q(iq,:,:),dyn_mat_qsym)
         dyn_mat_qsym(:,:) = dyn_mat_qsym(:,:)
         dyn_mat_q_full(iqfull,:,:) = dyn_mat_qsym

         ! Perform the actual FT onto the lattice vector grid
         call ft_dyn_direct(ft_lim,1,q_full,dyn_mat_qsym,fft_grid)
      end do 

      fft_grid(:,:,:)=fft_grid(:,:,:)/qpts%nkptf
      if (l_WSinterpol) then 
         ! unroll the FCM on supercell index
         ! for WS construction 
         iGrid = 1 
         do iz=ft_lim(1,3),ft_lim(2,3)
            do iy=ft_lim(1,2),ft_lim(2,2)
               do ix=ft_lim(1,1),ft_lim(2,1)
                  ! shift to storage indices (0-based)
                  nx = ix - ft_lim(1,1)
                  ny = iy - ft_lim(1,2)
                  nz = iz - ft_lim(1,3)
                  dyn_mat_r(nx,ny,nz,:,:)= fft_grid(iGrid,:,:)
                  iGrid = iGrid + 1 
               end do !ix
            end do !iy
         end do !iz
      end if 
   END SUBROUTINE

   SUBROUTINE ft_dyn_direct(ft_lim,isn,bqpt,dyn_mat_q,dyn_mat_r)
      INTEGER, INTENT(IN) :: ft_lim(2,3), isn
      REAL,    INTENT(IN) :: bqpt(3)

      COMPLEX,INTENT(INOUT) :: dyn_mat_q(:,:)
      COMPLEX,INTENT(INOUT) :: dyn_mat_r(:,:,:)

      INTEGER :: iGrid, ix, iy, iz
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

   subroutine ft_fcm_weight(ft_lim,bigBoxLim,weights,bqpt,dyn_mat_q,dyn_mat_r)
      ! Fourier transform for a fourier transform from a bigger box
      ! make use of weights to enforce correct periodicity 

      integer, intent(in) :: ft_lim(2,3), bigBoxlim(2,3)
      real,    intent(in) :: bqpt(3)
      real,    intent(in) :: weights(:)

      complex,intent(inout) :: dyn_mat_q(:,:)
      complex,intent(in) :: dyn_mat_r(0:,0:,0:,:,:)

      integer :: iGrid, ix, iy, iz, nx, ny, nz 
      real    :: phas
      complex :: phase_fac

      iGrid=0
      do iz=bigBoxlim(1,3),bigBoxlim(2,3)
         do iy=bigBoxlim(1,2),bigBoxlim(2,2)
            do ix=bigBoxlim(1,1),bigBoxlim(2,1)
               iGrid = iGrid+1
               phas=-1*tpi_const*(bqpt(1)*ix+bqpt(2)*iy+bqpt(3)*iz)
               phase_fac=cmplx(cos(phas),sin(phas))
               ! map to smaller box using ft_lim bounds
               nx = ft_lim(1,1) + modulo(ix - ft_lim(1,1), ft_lim(2,1) - ft_lim(1,1) + 1)
               ny = ft_lim(1,2) + modulo(iy - ft_lim(1,2), ft_lim(2,2) - ft_lim(1,2) + 1)
               nz = ft_lim(1,3) + modulo(iz - ft_lim(1,3), ft_lim(2,3) - ft_lim(1,3) + 1)
               ! shift to storage indices (0-based)
               nx = nx - ft_lim(1,1)
               ny = ny - ft_lim(1,2)
               nz = nz - ft_lim(1,3)
               dyn_mat_q(:,:)    = dyn_mat_q(:,:)    + phase_fac*weights(iGrid)*dyn_mat_r(nx,ny,nz,:,:)
            end do 
         end do 
      end do 
   end subroutine 

   SUBROUTINE rotate_dynmat(atoms,sym,isym,mrot,invmrot,l_inv,amat,bqpt,dyn,dyn_mat_qsym)
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
      !! $$f(\alpha,\beta)=exp(ix),x=\boldsymbol{q}_{\mathrm{red}}\cdot(\boldsymbol{\tau}_{\beta'}-\boldsymbol{\tau}_{\alpha'})-\boldsymbol{q}_{\mathrm{red}}\cdot(\boldsymbol{\tau}_{\beta}-\boldsymbol{\tau}_{\alpha})$$
      !! which can be written as a product
      !! $$f(\alpha,\beta) = f(\alpha)f^{*}(\beta), f(\alpha)=exp(i(\boldsymbol{q}_{\mathrm{red}}\cdot\boldsymbol{\tau}_{\alpha}-\boldsymbol{q}_{\mathrm{red}}\cdot\boldsymbol{\tau}_{\alpha'})).$$
      !! The real space rotation is related to the rotation matrix of the symmetry operation by the Bravais matrix of the system
      !! $$\underline{B}=\underline{A}\underline{S}^{-1}\underline{A}^{-1}$$
      USE m_inv3
      type(t_atoms), INTENT(IN)    :: atoms
      type(t_sym),   INTENT(IN)    :: sym
      INTEGER,       INTENT(IN)    :: isym
      INTEGER,          INTENT(IN)    :: mrot(3,3)
      INTEGER,          INTENT(IN)    :: invmrot(3,3) 
      LOGICAL,       INTENT(IN)    :: l_inv
      REAL,          INTENT(IN)    :: amat(3,3)
      REAL,          INTENT(IN)    :: bqpt(3)
      COMPLEX,       INTENT(IN)    :: dyn(:,:)
      COMPLEX,       INTENT(INOUT) :: dyn_mat_qsym(:,:)
      
      INTEGER :: iAtom, jAtom , iAtom_map, jAtom_map
      INTEGER :: iAlpha, iBeta, iAlpha_map, iBeta_map
      REAL    :: invamat(3,3), phas, det , pos_map(3)
      COMPLEX :: brot(3,3), temp_mat_1(3,3), temp_mat_2(3,3)
      COMPLEX :: phase_fac

      INTEGER :: map(atoms%nat)
      COMPLEX :: phase_map(atoms%nat)

      CALL inv3(amat,invamat,det)
      temp_mat_1 = MATMUL(invmrot,invamat)
      brot = MATMUL(amat,temp_mat_1)
      
      ! Find the q vector in the full BZ
      !q_full = MATMUL(transpose(mrot),bqpt)
      
      ! Calculate the array of phases
      DO iAtom = 1, atoms%nat
         jAtom = sym%mapped_atom(isym,iAtom)
         ! Create list to which atom we map 
         map(iAtom)=jAtom
         ! Find rotated atom position
         pos_map = matmul(mrot, atoms%taual(:,iAtom))
         IF (l_inv) pos_map = -1 * pos_map ! Inversion does not exist in real space, 
                                           ! we have to introduce a minus that negates the minus of mrot
         ! Calculate the phase factor f(alpha)
         phas= tpi_const*(dot_product(bqpt(:),atoms%taual(:,jAtom)) - dot_product(bqpt(:), pos_map(:))) 
         phase_fac=cmplx(cos(phas),sin(phas))
         phase_map(iAtom)=phase_fac
      END DO
 
      IF (l_inv) phase_map = conjg(phase_map) ! inversion maps q -> -q which results in a complex conjugate of the phase
      ! Transform the dynamical matrix from the representative atom and q vector to the unfolded ones
      DO iAtom=1, atoms%nat
         iAlpha = 3*(iAtom-1)
         iAlpha_map = 3*(map(iAtom)-1)
         iAtom_map = map(iAtom)
         DO jAtom=1, atoms%nat
            iBeta = 3*(jAtom-1)
            iBeta_map = 3*(map(jAtom)-1)
            jAtom_map = map(jAtom)
            temp_mat_1 = dyn(iAlpha+1:iAlpha+3,iBeta+1:iBeta+3)
            if (l_inv) temp_mat_1 = conjg(temp_mat_1) ! inversion maps q -> -q  , which results into a conjugate
            temp_mat_2 = MATMUL(brot,temp_mat_1)
            temp_mat_1 = MATMUL(temp_mat_2,TRANSPOSE(brot))
            phase_fac=phase_map(iAtom_map)*conjg(phase_map(jAtom_map))

            dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) &
          = dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) + phase_fac*temp_mat_1
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE ift_dyn(atoms,qpts,ft_lim,bigBox_lim,weights,bqpt,dyn_mat_r,fft_grid,dyn_mat_q,l_WSinterpol)
      !! Transforms the dynamical matrix on a real space lattice vector grid
      !! (--> mass-normalized FCM) back onto a specific q vector provided as
      !! input (bqpt) by the inverse Fourier Transformation as compared to
      !! SUBROUTINE ft_dyn.
      type(t_atoms), intent(in) :: atoms
      type(t_kpts),  intent(in) :: qpts
      integer,  intent(in) :: ft_lim(2,3), bigBox_lim(2,3)
      real,     intent(in) :: weights(:) 
      real,          intent(in) :: bqpt(3)
      complex,       intent(in) :: dyn_mat_r(0:,0:,0:,:,:)
      complex,       intent(inout) :: fft_grid(:,:,:)
      complex, allocatable, intent(out) :: dyn_mat_q(:,:)
      logical, intent(in) :: l_WSinterpol

      integer :: dyn_dim 

      dyn_dim = 3*atoms%nat

      allocate(dyn_mat_q(dyn_dim,dyn_dim))
      dyn_mat_q(:,:) = cmplx(0.0,0.0)

      if (l_WSinterpol) then 
         call ft_fcm_weight(ft_lim,bigBox_lim,weights,bqpt,dyn_mat_q,dyn_mat_r)
      else 
         CALL ft_dyn_direct(ft_lim,-1,bqpt,dyn_mat_q,fft_grid)
      end if 

   END SUBROUTINE

   SUBROUTINE make_sym_list(sym, bqpt, sym_count, sym_list)
      TYPE(t_sym), INTENT(IN)    :: sym
      REAL,        INTENT(IN)    :: bqpt(3)
      INTEGER,     INTENT(OUT)   :: sym_count
      INTEGER,     INTENT(INOUT) :: sym_list(:)
      
      INTEGER :: iSym
      
      sym_count = 0
      sym_list = 0
      DO iSym = 1, sym%nop
         IF (norm2(bqpt-MATMUL(bqpt,sym%mrot(:,:,iSym)))<1e-8) THEN
            sym_count = sym_count + 1
            sym_list(sym_count) = iSym
         END IF
      END DO
   END SUBROUTINE

   SUBROUTINE make_sym_dynvec(atoms, sym, amat, bqpt, iDtype, iDir, sym_count, sym_list, dynvec, sym_dynvec)
      USE m_inv3

      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_sym),   INTENT(IN)    :: sym
      REAL,          INTENT(IN)    :: amat(3,3), bqpt(3)
      INTEGER,       INTENT(IN)    :: iDtype, iDir, sym_count
      INTEGER,       INTENT(IN)    :: sym_list(sym_count)
      COMPLEX,       INTENT(IN)    :: dynvec(:)
      COMPLEX,       INTENT(INOUT) :: sym_dynvec(:,:,:)
     
      INTEGER :: iSym, iAtom, iRow, iCol
      REAL    :: phas, det, mrot(3,3), invmrot(3,3), invamat(3,3)
      COMPLEX :: brot(3,3), temp_mat_1(3,3), temp_mat_2(3,3)
      COMPLEX :: phase_fac, rotvec(3)

      iRow = 3 * (iDtype-1) + iDir

      DO iSym = 1, sym_count
         mrot = sym%mrot(:,:,sym_list(iSym))
         invmrot = sym%mrot(:,:,sym%invtab(sym_list(iSym)))
         CALL inv3(amat,invamat,det)
         temp_mat_1 = MATMUL(invmrot,invamat)
         brot = MATMUL(amat,temp_mat_1)
         DO iAtom = 1, atoms%nat
            iCol = 3 * (iAtom-1)

            phas = -tpi_const*(dot_product(bqpt(:),atoms%taual(:,iAtom)-atoms%taual(:,iDtype)))
            phase_fac = cmplx(cos(phas),sin(phas))

            rotvec = MATMUL(brot,dynvec(iCol+1:iCol+3))
            sym_dynvec(iCol+1:iCol+3,iRow,iSym) = phase_fac * rotvec
         END DO
      END DO 
   END SUBROUTINE

   SUBROUTINE cheat_dynmat(atoms, sym, amat, bqpt, iBetaPr, jDirPr, sym_count, sym_list, sym_dynvec, dynmat, sym_dynmat, l_cheated)
      USE m_inv3

      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_sym),   INTENT(IN)    :: sym
      REAL,          INTENT(IN)    :: amat(3,3), bqpt(3)
      INTEGER,       INTENT(IN)    :: iBetaPr, jDirPr, sym_count
      INTEGER,       INTENT(IN)    :: sym_list(sym_count)
      COMPLEX,       INTENT(IN)    :: sym_dynvec(:,:,:), dynmat(:,:)
      COMPLEX,       INTENT(INOUT) :: sym_dynmat(:,:)
      LOGICAL,       INTENT(OUT)   :: l_cheated

      INTEGER :: iSym, iDatom, iAtom, iRowPr, iColPr, iRow, iCol, iDir, jDir, iDirPr, kDir, iBeta, iAlpha, iAlphaPr, symsumcount, aAlpha(atoms%nat), iDone, iNeed, kNeed
      REAL    :: phas, det, mrot(3,3), invmrot(3,3), invamat(3,3)
      COMPLEX :: brot(3,3), temp_mat_1(3,3), temp_mat_2(3,3), symsum(3,3), symvec(3)
      COMPLEX :: phase_fac, rotvec(3), rhs_vecs(3,3 * (iBetaPr-1) + jDirPr - 1), rhs_vec_full(3)
      LOGICAL :: l_mapped, l_done(3), l_need(3)

      l_cheated = .FALSE.
      iRowPr = 3 * (iBetaPr-1) + jDirPr
      alphaPrloop: DO iAlphaPr = 1, atoms%nat
         iColPr = 3 * (iAlphaPr-1)
         DO iBeta = 1, iBetaPr
            iRow = 3 * (iBeta-1)
            alphaloop: DO iAlpha = 1, atoms%nat
               iCol = 3 * (iAlpha-1)
               symloop: DO iSym = 1, sym_count
                  IF (.NOT.(iBetaPr==sym%mapped_atom(sym_list(iSym),iBeta))) CYCLE
                  IF (.NOT.(iAlphaPr==sym%mapped_atom(sym_list(iSym),iAlpha))) CYCLE
                  ! Get all rhs vectors that are possible
                  rhs_vecs = sym_dynvec(iCol+1:iCol+3,:iRowPr-1,iSym)
                  phas = tpi_const*(dot_product(bqpt(:),atoms%taual(:,iAlphaPr)-atoms%taual(:,iBetaPr)))
                  phase_fac = cmplx(cos(phas),sin(phas))
                  
                  invmrot = sym%mrot(:,:,sym%invtab(sym_list(iSym)))
                  CALL inv3(amat,invamat,det)
                  temp_mat_1 = MATMUL(invmrot,invamat)
                  brot = MATMUL(amat,temp_mat_1)
                  
                  DO jDir = 1, 3
                     IF (iRow+jDir>iRowPr) CYCLE symloop
                     rhs_vec_full = phase_fac*rhs_vecs(:, iRow + jDir)
                     !write(*,*) "---------------"
                     !write(*,*) rhs_vec_full
                     symvec = brot(:, jDir)
                     !write(*,*) symvec
                  
                     l_done = .FALSE.
                     iDone = 0
                     l_need = .TRUE.
                     iNeed = 3
                     kNeed = 0
                     DO kDir = 1, 3
                        !IF (ABS(symvec(kDir))>1e-8.AND.ANY((ABS(dynmat(3 * (iBetaPr-1) + kDir,:))<1e-15)))  
                        IF (ALL(ABS(dynmat(3 * (iBetaPr-1) + kDir,iColPr+1:iColPr+3))>1e-15)) THEN 
                           l_done(kDir) = .TRUE.
                           iDone = iDone + 1
                           l_need(kDir) = .FALSE.
                           iNeed = iNeed - 1
                           CYCLE
                        END IF
                        IF (ABS(symvec(kDir))<1e-8) THEN
                           l_need(kDir) = .FALSE.
                            iNeed = iNeed - 1
                            CYCLE
                        END IF
                        kNeed = kDir
                     END DO
                     !write(*,*) iDone, l_done
                     !write(*,*) iNeed, l_need, kNeed
                     IF (iDone==0.OR.iDone==3) CYCLE
                     IF (iNeed/=1) CYCLE
                     DO kDir = 1, 3
                        IF (.NOT.l_need(kDir)) rhs_vec_full = rhs_vec_full - symvec(kDir)*dynmat(3 * (iBetaPr-1) + kDir,iColPr+1:iColPr+3)
                     END DO
                     IF (SQRT(REAL(DOT_PRODUCT(rhs_vec_full,rhs_vec_full)))<1e-15) CYCLE
                     !write(*,*) "Newrhs:", rhs_vec_full
                     !write(*,*) "Thesym:", symvec(kNeed)
                     IF (ABS(symvec(kNeed))>1e-8) THEN 
                        sym_dynmat(iRowPr-jDirPr+kNeed,iColPr+1:iColPr+3) = rhs_vec_full/symvec(kNeed)
                        l_cheated = .TRUE.
                     ELSE
                        CYCLE
                     END IF
                     !write(*,*) "Cheat: ", sym_dynmat(iRowPr-jDirPr+kNeed,:)
                     IF (ALL(ABS(sym_dynmat(iRowPr-jDirPr+kNeed,iColPr+1:iColPr+3))>1e-15)) CYCLE alphaPrloop
                  END DO
               END DO symloop
            END DO alphaloop
         END DO
      END DO alphaPrloop
   END SUBROUTINE

END MODULE m_dfpt_dynmat_sym
