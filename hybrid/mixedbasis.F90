!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates the mixed basis set used to evaluate the  !
! exchange term in HF/hybrid functional calculations or EXX           !
! calculations. In the latter case a second mixed basis set is setup  !
! for the OEP integral equation.                                      !
! In all cases the mixed basis consists of IR plane waves             !
!                                                                     !
! IR:                                                                 !
!    M_{\vec{k},\vec{G}} = 1/\sqrt{V} \exp{i(\vec{k}+\vec{G})}        !
!                                                                     !
! which are zero in the MT spheres and radial functions times         !
! spherical harmonics in the MT spheres                               !
!                                                                     !
! MT:                                                                 !
!     a                a                                              !
!    M              = U   Y                                           !
!     PL               PL  LM                                         !
!                                                                     !
!            where     a    a  a                                      !
!                     U  = u  u                                       !
!                      PL   pl  p'l'                                  !
!                                                                     !
!               and    L \in {|l-l'|,...,l+l'}                        !
!                                                                     !
!               and    P counts the different combinations of         !
!                      pl,p'l' which contribute to L                  !
!                                                                     !
!                                               M.Betzinger (09/07)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_mixedbasis

CONTAINS

   SUBROUTINE mixedbasis(atoms, kpts, input, cell, xcpot, mpbasis, hybrid, enpara, mpi, v)

      USE m_judft
      USE m_loddop, ONLY: loddop
      USE m_intgrf, ONLY: intgrf_init, intgrf
      use m_rorder, only: rorderpf
      USE m_hybrid_core
      USE m_wrapper
      USE m_eig66_io
      USE m_types

      IMPLICIT NONE

      TYPE(t_xcpot_inbuild), INTENT(IN)    :: xcpot
      TYPE(t_mpi), INTENT(IN)    :: mpi
      TYPE(t_mpbasis), intent(inout)  :: mpbasis
      TYPE(t_hybrid), INTENT(INOUT) :: hybrid
      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_potden), INTENT(IN)    :: v


      ! local type variables
      TYPE(t_usdus)                   ::  usdus

      ! local scalars
      INTEGER                         ::  jspin, itype, l1, l2, l, n_radbasfn, n1, n2, nn
      INTEGER                         ::  m, nk, i_basfn, i, j, n_grid_pt
      REAL                            ::  rdum, rdum1, norm

      ! - local arrays -
      INTEGER, ALLOCATABLE             ::  ihelp(:)

      REAL                            ::  bashlp(atoms%jmtd)


      REAL, ALLOCATABLE               ::  olap(:,:), work(:), eig(:), eigv(:,:)
      REAL, ALLOCATABLE               ::  bas1(:,:,:,:,:), bas2(:,:,:,:,:)
      REAL, ALLOCATABLE               ::  basmhlp(:,:,:,:)
      REAL, ALLOCATABLE               ::  gridf(:,:), vr0(:,:,:)

      LOGICAL, ALLOCATABLE            ::  selecmat(:,:,:,:)
      LOGICAL, ALLOCATABLE            ::  seleco(:,:), selecu(:,:)

      CHARACTER, PARAMETER            :: lchar(0:38) = (/'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', &
                                                         'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', &
                                                         'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x'/)


      IF (mpi%irank == 0) WRITE (6, '(//A,I2,A)') '### subroutine: mixedbasis ###'

      IF (xcpot%is_name("exx")) CALL judft_error("EXX is not implemented in this version", calledby='mixedbasis')

      ! Deallocate arrays which might have been allocated in a previous run of this subroutine
      IF (ALLOCATED(mpbasis%ngptm)) deallocate(mpbasis%ngptm)
      IF (ALLOCATED(mpbasis%num_rad_bas_fun)) deallocate(mpbasis%num_rad_bas_fun)
      IF (ALLOCATED(mpbasis%gptm_ptr)) deallocate(mpbasis%gptm_ptr)
      IF (ALLOCATED(mpbasis%gptm)) deallocate(mpbasis%gptm)
      IF (ALLOCATED(hybrid%radbasfn_mt)) deallocate(hybrid%radbasfn_mt)

      CALL usdus%init(atoms, input%jspins)
      call hybrid%set_num_radfun_per_l(atoms)

      ! initialize gridf for radial integration
      CALL intgrf_init(atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, gridf)

      allocate(vr0(atoms%jmtd, atoms%ntype, input%jspins), source=0.0)

      vr0(:,:,:) = v%mt(:,0, :,:)

      ! calculate radial basisfunctions u and u' with
      ! the spherical part of the potential vr0 and store them in
      ! bas1 = large component ,bas2 = small component

      call gen_bas_fun(atoms, enpara, gridf, input, hybrid, mpi, vr0, usdus, bas1, bas2)

      ! - - - - - - SETUP OF THE MIXED BASIS IN THE IR - - - - - - -

      ! construct G-vectors with cutoff smaller than gcutm
      call mpbasis%gen_gvec(cell, kpts, mpi)

      ! - - - - - - - - Set up MT product basis for the non-local exchange potential  - - - - - - - - - -

      IF (mpi%irank == 0) THEN
         WRITE (6, '(A)') 'MT product basis for non-local exchange potential:'
         WRITE (6, '(A)') 'Reduction due to overlap (quality of orthonormality, should be < 1.0E-06)'
      END IF

      allocate(mpbasis%num_rad_bas_fun(0:maxval(hybrid%lcutm1), atoms%ntype))
      allocate(seleco(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd))
      allocate(selecu(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd))
      mpbasis%num_rad_bas_fun = 0    !!! 01/12/10 jij%M.b.

      ! determine maximal indices of (radial) mixed-basis functions (->num_rad_bas_fun)
      ! (will be reduced later-on due to overlap)
      hybrid%max_indx_p_1 = 0
      DO itype = 1, atoms%ntype
         seleco = .FALSE.
         selecu = .FALSE.
         seleco(1, 0:hybrid%select1(1, itype)) = .TRUE.
         selecu(1, 0:hybrid%select1(3, itype)) = .TRUE.
         seleco(2, 0:hybrid%select1(2, itype)) = .TRUE.
         selecu(2, 0:hybrid%select1(4, itype)) = .TRUE.

         ! include local orbitals
         IF (maxval(hybrid%num_radfun_per_l) >= 3) THEN
            seleco(3:,:) = .TRUE.
            selecu(3:,:) = .TRUE.
         END IF

         DO l = 0, hybrid%lcutm1(itype)
            n_radbasfn = 0
            M = 0

            !
            ! valence * valence
            !
            if(.not. allocated(selecmat)) then
               allocate(selecmat(maxval(hybrid%num_radfun_per_l), &
                                 0:atoms%lmaxd, &
                                 maxval(hybrid%num_radfun_per_l), &
                                 0:atoms%lmaxd))
            endif
            selecmat = calc_selecmat(atoms, hybrid, seleco, selecu)

            DO l1 = 0, atoms%lmax(itype)
               DO l2 = 0, atoms%lmax(itype)
                  IF (l >= ABS(l1 - l2) .AND. l <= l1 + l2) THEN
                     DO n1 = 1, hybrid%num_radfun_per_l(l1, itype)
                        DO n2 = 1, hybrid%num_radfun_per_l(l2, itype)
                           M = M + 1
                           IF (selecmat(n1, l1, n2, l2)) THEN
                              n_radbasfn = n_radbasfn + 1
                              selecmat(n2, l2, n1, l1) = .FALSE. ! prevent double counting of products (a*b = b*a)
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
            IF (n_radbasfn == 0 .AND. mpi%irank == 0) &
               WRITE (6, '(A)') 'mixedbasis: Warning!  No basis-function product of '//lchar(l)// &
               '-angular momentum defined.'
            hybrid%max_indx_p_1 = MAX(hybrid%max_indx_p_1, M)
            mpbasis%num_rad_bas_fun(l, itype) = n_radbasfn*input%jspins
         END DO
      END DO

      allocate(hybrid%radbasfn_mt(atoms%jmtd,&
                            maxval(mpbasis%num_rad_bas_fun), &
                            0:maxval(hybrid%lcutm1), &
                            atoms%ntype), source=0.0)

      ! Define product bases and reduce them according to overlap

      DO itype = 1, atoms%ntype
         seleco = .FALSE.
         selecu = .FALSE.
         seleco(1, 0:hybrid%select1(1, itype)) = .TRUE.
         selecu(1, 0:hybrid%select1(3, itype)) = .TRUE.
         seleco(2, 0:hybrid%select1(2, itype)) = .TRUE.
         selecu(2, 0:hybrid%select1(4, itype)) = .TRUE.
         ! include lo's
         IF (maxval(hybrid%num_radfun_per_l) >= 3) THEN
            seleco(3:,:) = .TRUE.
            selecu(3:,:) = .TRUE.
         END IF

         n_grid_pt = atoms%jri(itype)
         DO l = 0, hybrid%lcutm1(itype)
            n_radbasfn = mpbasis%num_rad_bas_fun(l, itype)
            ! allow for zero product-basis functions for
            ! current l-quantum number
            IF (n_radbasfn == 0) THEN
               IF (mpi%irank == 0) WRITE (6, '(6X,A,'':   0 ->   0'')') lchar(l)
               CYCLE
            END IF

            ! set up the overlap matrix
            allocate(olap(n_radbasfn, n_radbasfn), source=0.0)
            allocate(eigv(n_radbasfn, n_radbasfn), source=0.0)
            allocate(work(3*n_radbasfn), source=0.0)
            allocate(eig(n_radbasfn), source=0.0)
            allocate(ihelp(n_radbasfn), source=1)
            i_basfn = 0

            ! valence*valence
            selecmat =  calc_selecmat(atoms, hybrid, seleco, selecu)

            DO l1 = 0, atoms%lmax(itype)
               DO l2 = 0, atoms%lmax(itype)
                  IF (l >= ABS(l1 - l2) .AND. l <= l1 + l2) THEN
                     DO n1 = 1, hybrid%num_radfun_per_l(l1, itype)
                        DO n2 = 1, hybrid%num_radfun_per_l(l2, itype)

                           IF (selecmat(n1, l1, n2, l2)) THEN
                              DO jspin = 1, input%jspins
                                 i_basfn = i_basfn + 1
                                 IF (i_basfn > n_radbasfn) call judft_error('got too many product functions', hint='This is a BUG, please report', calledby='mixedbasis')

                                 hybrid%radbasfn_mt(:n_grid_pt, i_basfn, l, itype) &
                                    = (   bas1(:n_grid_pt, n1, l1, itype, jspin) &
                                        * bas1(:n_grid_pt, n2, l2, itype, jspin) &
                                        + bas2(:n_grid_pt, n1, l1, itype, jspin) &
                                        * bas2(:n_grid_pt, n2, l2, itype, jspin) &
                                      ) / atoms%rmsh(:n_grid_pt, itype)

                                 !normalize radbasfn_mt
                                 norm = calc_radbas_norm(atoms, hybrid, i_basfn, itype, gridf)
                                 hybrid%radbasfn_mt(:n_grid_pt, i_basfn, l, itype) &
                                 = hybrid%radbasfn_mt(:n_grid_pt, i_basfn, l, itype)/norm

                              END DO !jspin
                              ! prevent double counting of products (a*b = b*a)
                              selecmat(n2, l2, n1, l1) = .FALSE.
                           END IF
                        END DO !n2
                     END DO !n1
                  ENDIF
               END DO !l2
            END DO  !l1

            IF (i_basfn /= n_radbasfn) call judft_error('counting error for product functions', hint='This is a BUG, please report', calledby='mixedbasis')

            ! In order to get rid of the linear dependencies in the
            ! radial functions radbasfn_mt belonging to fixed l and itype
            ! the overlap matrix is diagonalized and those eigenvectors
            ! with a eigenvalue greater then hybrid%tolerance1 are retained

            ! Calculate overlap
            DO n2 = 1, n_radbasfn
               DO n1 = 1, n2
                  olap(n1, n2) = intgrf(hybrid%radbasfn_mt(:,n1, l, itype)*hybrid%radbasfn_mt(:,n2, l, itype), &
                                        atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
                  olap(n2, n1) = olap(n1, n2)
               END DO
            END DO

            ! Diagonalize
            CALL diagonalize(eigv, eig, olap)

            ! Get rid of linear dependencies (eigenvalue <= hybrid%tolerance1)
            nn = 0
            DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
               IF (eig(i) > hybrid%tolerance1) THEN
                  nn = nn + 1
                  ihelp(nn) = i
               END IF
            END DO
            mpbasis%num_rad_bas_fun(l, itype) = nn
            eig = eig(ihelp)
            eigv(:,:) = eigv(:,ihelp)

            DO i = 1, n_grid_pt
               hybrid%radbasfn_mt(i, 1:nn, l, itype) = MATMUL(hybrid%radbasfn_mt(i, 1:n_radbasfn, l, itype), eigv(:,1:nn))/SQRT(eig(:nn))
            END DO

            ! Add constant function to l=0 basis and then do a Gram-Schmidt orthonormalization
            IF (l == 0) THEN

               ! Check if radbasfn_mt must be reallocated
               IF (nn + 1 > SIZE(hybrid%radbasfn_mt, 2)) THEN
                  allocate(basmhlp(atoms%jmtd, nn + 1, 0:maxval(hybrid%lcutm1), atoms%ntype))
                  basmhlp(:,1:nn, :,:) = hybrid%radbasfn_mt
                  deallocate(hybrid%radbasfn_mt)
                  allocate(hybrid%radbasfn_mt(atoms%jmtd, nn + 1, 0:maxval(hybrid%lcutm1), atoms%ntype))
                  hybrid%radbasfn_mt(:,1:nn, :,:) = basmhlp(:,1:nn, :,:)
                  deallocate(basmhlp)
               END IF

               hybrid%radbasfn_mt(:n_grid_pt, nn + 1, 0, itype) = atoms%rmsh(:n_grid_pt, itype)/SQRT(atoms%rmsh(n_grid_pt, itype)**3/3)
               DO i = nn, 1, -1
                  DO j = i + 1, nn + 1
                     hybrid%radbasfn_mt(:n_grid_pt, i, 0, itype) = hybrid%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                                                      - intgrf(hybrid%radbasfn_mt(:n_grid_pt, i, 0, itype)*hybrid%radbasfn_mt(:n_grid_pt, j, 0, itype), &
                                                               atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf) &
                                                      *hybrid%radbasfn_mt(:n_grid_pt, j, 0, itype)

                  END DO
                  hybrid%radbasfn_mt(:n_grid_pt, i, 0, itype) = hybrid%radbasfn_mt(:n_grid_pt, i, 0, itype) &
                                                   /SQRT(intgrf(hybrid%radbasfn_mt(:n_grid_pt, i, 0, itype)**2, &
                                                                atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf))
               END DO
               nn = nn + 1
               deallocate(olap)
               allocate(olap(nn, nn))
               mpbasis%num_rad_bas_fun(l, itype) = nn
            END IF

            ! Check orthonormality of product basis
            rdum = 0
            DO i = 1, nn
               DO j = 1, i
                  rdum1 = intgrf(hybrid%radbasfn_mt(:,i, l, itype)*hybrid%radbasfn_mt(:,j, l, itype), &
                                 atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

                  IF (i == j) THEN
                     rdum1 = ABS(1 - rdum1)
                     rdum = rdum + rdum1**2
                  ELSE
                     rdum1 = ABS(rdum1)
                     rdum = rdum + 2*rdum1**2
                  END IF

                  IF (rdum1 > 1e-6) THEN
                     IF (mpi%irank == 0) THEN
                        WRITE (6, '(A)') 'mixedbasis: Bad orthonormality of ' &
                           //lchar(l)//'-product basis. Increase tolerance.'
                        WRITE (6, '(12X,A,F9.6,A,2(I3.3,A))') 'Deviation of', &
                           rdum1, ' occurred for (', i, ',', j, ')-overlap.'
                     END IF
                     CALL judft_error("Bad orthonormality of product basis", hint='Increase tolerance', calledby='mixedbasis')
                  END IF

               END DO
            END DO

            IF (mpi%irank == 0) THEN
               WRITE (6, '(6X,A,I4,'' ->'',I4,''   ('',ES8.1,'' )'')') lchar(l)//':', n_radbasfn, nn, SQRT(rdum)/nn
            END IF

            deallocate(olap, eigv, work, eig, ihelp)

         END DO !l
         IF (mpi%irank == 0) WRITE (6, '(6X,A,I7)') 'Total:', SUM(mpbasis%num_rad_bas_fun(0:hybrid%lcutm1(itype), itype))
      END DO ! itype

      allocate(basmhlp(atoms%jmtd, maxval(mpbasis%num_rad_bas_fun), 0:maxval(hybrid%lcutm1), atoms%ntype))
      basmhlp(1:atoms%jmtd, 1:maxval(mpbasis%num_rad_bas_fun), 0:maxval(hybrid%lcutm1), 1:atoms%ntype) &
         = hybrid%radbasfn_mt(1:atoms%jmtd, 1:maxval(mpbasis%num_rad_bas_fun), 0:maxval(hybrid%lcutm1), 1:atoms%ntype)
      deallocate(hybrid%radbasfn_mt)
      allocate(hybrid%radbasfn_mt(atoms%jmtd, maxval(mpbasis%num_rad_bas_fun), 0:maxval(hybrid%lcutm1), atoms%ntype))
      hybrid%radbasfn_mt = basmhlp

      deallocate(basmhlp, seleco, selecu, selecmat)

      !
      ! now we build linear combinations of the radial functions
      ! such that they possess no moment except one radial function in each l-channel
      !
      IF (mpi%irank == 0) THEN
         WRITE (6, '(/,A,/,A)') 'Build linear combinations of radial '// &
            'functions in each l-channel,', &
            'such that they possess no multipolmoment'// &
            ' except the last function:'

         WRITE (6, '(/,17x,A)') 'moment  (quality of orthonormality)'
      END IF
      DO itype = 1, atoms%ntype
         n_grid_pt = atoms%jri(itype)

         IF (atoms%ntype > 1 .AND. mpi%irank == 0) WRITE (6, '(6X,A,I3)') 'Atom type', itype

         DO l = 0, hybrid%lcutm1(itype)
            ! determine radial function with the largest moment
            ! this function is used to build the linear combinations
            rdum = 0
            DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
               rdum1 = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*hybrid%radbasfn_mt(:n_grid_pt, i, l, itype), &
                              atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
               IF (ABS(rdum1) > rdum) THEN
                  n_radbasfn = i
                  rdum = rdum1
               END IF
            END DO

            ! rearrange order of radial functions such that the last function possesses the largest moment
            j = 0
            bashlp(:n_grid_pt) = hybrid%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype)
            DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
               IF (i == n_radbasfn) CYCLE
               j = j + 1
               hybrid%radbasfn_mt(:n_grid_pt, j, l, itype) = hybrid%radbasfn_mt(:n_grid_pt, i, l, itype)
            END DO
            hybrid%radbasfn_mt(:n_grid_pt, mpbasis%num_rad_bas_fun(l, itype), l, itype) = bashlp(:n_grid_pt)

         END DO

         DO l = 0, hybrid%lcutm1(itype)
            IF (mpi%irank == 0) WRITE (6, '(6X,A)') lchar(l)//':'

            IF (mpbasis%num_rad_bas_fun(l, itype) == 0) THEN
               IF (mpi%irank == 0) WRITE (6, '(6X,A,'':   0 ->    '')') lchar(l)
               CYCLE
            END IF

            n_radbasfn = mpbasis%num_rad_bas_fun(l, itype)
            DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
               IF (i == n_radbasfn) CYCLE
               ! calculate moment of radial function i
               rdum1 = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*hybrid%radbasfn_mt(:n_grid_pt, i, l, itype), &
                              atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

               rdum = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*hybrid%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype), &
                             atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

               bashlp(:n_grid_pt) = hybrid%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype)

               IF (SQRT(rdum**2 + rdum1**2) <= 1E-06 .AND. mpi%irank == 0) &
                  WRITE (6, *) 'Warning: Norm is smaller thann 1E-06!'

               ! change function n_radbasfn such that n_radbasfn is orthogonal to i
               ! since the functions radbasfn_mt have been orthogonal on input
               ! the linear combination does not destroy the orthogonality to the residual functions
               hybrid%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype) = rdum/SQRT(rdum**2 + rdum1**2)*bashlp(:n_grid_pt) &
                                                + rdum1/SQRT(rdum**2 + rdum1**2)*hybrid%radbasfn_mt(:n_grid_pt, i, l, itype)

               ! combine basis function i and n_radbasfn so that they possess no momemt
               hybrid%radbasfn_mt(:n_grid_pt, i, l, itype) = rdum1/SQRT(rdum**2 + rdum1**2)*bashlp(:n_grid_pt) &
                                                - rdum/SQRT(rdum**2 + rdum1**2)*hybrid%radbasfn_mt(:n_grid_pt, i, l, itype)

               rdum1 = intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*hybrid%radbasfn_mt(:n_grid_pt, i, l, itype), &
                              atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

               IF (rdum1 > 1E-10) call judft_error('moment of radial function does not vanish', calledby='mixedbasis')

               IF (mpi%irank == 0) WRITE (6, '(6x,I4,'' ->  '',ES8.1)') i, rdum1
            END DO

            ! test orthogonality
            rdum = 0
            DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
               DO j = 1, mpbasis%num_rad_bas_fun(l, itype)
                  rdum1 = intgrf(hybrid%radbasfn_mt(:n_grid_pt, i, l, itype)*hybrid%radbasfn_mt(:n_grid_pt, j, l, itype), &
                                 atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
                  IF (i /= j) THEN
                     rdum = rdum + rdum1
                  ELSE
                     rdum = rdum + ABS(1 - rdum1)
                  END IF
               END DO
            END DO
            IF (mpi%irank == 0) &
               WRITE (6, '(6x,I4,'' ->'',f10.5,''   ('',ES8.1,'' )'')') n_radbasfn, &
               intgrf(atoms%rmsh(:n_grid_pt, itype)**(l + 1)*hybrid%radbasfn_mt(:n_grid_pt, n_radbasfn, l, itype), atoms%jri, &
                      atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf), rdum
         END DO
      END DO

      DO itype = 1, atoms%ntype
         IF (ANY(mpbasis%num_rad_bas_fun(0:hybrid%lcutm1(itype), itype) == 0)) call judft_error('any mpbasis%num_rad_bas_fun eq 0', calledby='mixedbasis')
      END DO

      !count basis functions
      hybrid%nbasp = 0
      DO itype = 1, atoms%ntype
         DO i = 1, atoms%neq(itype)
            DO l = 0, hybrid%lcutm1(itype)
               DO M = -l, l
                  DO j = 1, mpbasis%num_rad_bas_fun(l, itype)
                     hybrid%nbasp = hybrid%nbasp + 1
                  END DO
               END DO
            END DO
         END DO
      END DO
      hybrid%maxbasm1 = hybrid%nbasp + maxval(mpbasis%ngptm)
      DO nk = 1, kpts%nkptf
         hybrid%nbasm(nk) = hybrid%nbasp + mpbasis%ngptm(nk)
      END DO

      hybrid%maxlmindx = MAXVAL([(SUM([(hybrid%num_radfun_per_l(l, itype)*(2*l + 1), l=0, atoms%lmax(itype))]), itype=1, atoms%ntype)])

   END SUBROUTINE mixedbasis

   subroutine gen_bas_fun(atoms, enpara, gridf, input, hybrid, mpi, vr0, usdus, bas1, bas2)
      use m_judft
      use m_types
      USE m_radfun, ONLY: radfun
      USE m_radflo, ONLY: radflo
      USE m_intgrf,   ONLY: intgrf
      implicit NONE
      type(t_atoms), intent(in)        :: atoms
      type(t_enpara), intent(in)       :: enpara
      type(t_input), intent(in)        :: input
      type(t_hybrid), intent(in)       :: hybrid
      type(t_mpi), intent(in)          :: mpi
      type(t_usdus), intent(inout)     :: usdus

      REAL, ALLOCATABLE, INTENT(INOUT) :: bas1(:,:,:,:,:), bas2(:,:,:,:,:)
      REAL, intent(in)                 :: vr0(:,:,:), gridf(:,:)

      REAL    ::   u(atoms%jmtd, 2, 0:atoms%lmaxd)
      REAL    ::  du(atoms%jmtd, 2, 0:atoms%lmaxd)
      REAL    :: flo(atoms%jmtd, 2, atoms%nlod)

      REAL    :: uuilon(atoms%nlod, atoms%ntype)
      REAL    :: duilon(atoms%nlod, atoms%ntype)
      REAL    :: ulouilopn(atoms%nlod, atoms%nlod, atoms%ntype)
      REAL    :: wronk, norm

      INTEGER :: itype, jspin, i, l, ilo, ok
      INTEGER :: n_grid_pt, noded, nodem
      u  = 0.0
      du = 0.0

      ! this is 5-D array. it could cause Problems in bigger systems
      allocate(bas1(atoms%jmtd,    &
                    maxval(hybrid%num_radfun_per_l), &
                    0:atoms%lmaxd, &
                    atoms%ntype,   &
                    input%jspins),   source=0.0, stat=ok)
      if(ok /= 0) call judft_error("Can't allocate bas1 array. Stat= " // int2str(ok))

      allocate(bas2, source=bas1, stat=ok)
      if(ok /= 0) call judft_error("Can't allocate bas1 array. Stat= " // int2str(ok))

      DO itype = 1, atoms%ntype
         n_grid_pt = atoms%jri(itype) ! number of radial gridpoints
         DO jspin = 1, input%jspins
            DO l = 0, atoms%lmax(itype)
               CALL radfun(l, itype, jspin, enpara%el0(l, itype, jspin), vr0(:,itype, jspin), atoms, &
                           u(:,:,l), du(:,:,l), usdus, nodem, noded, wronk)
            END DO
            bas1(1:n_grid_pt, 1, 0:atoms%lmaxd, itype, jspin)  = u(1:n_grid_pt, 1, 0:atoms%lmaxd)
            bas2(1:n_grid_pt, 1, 0:atoms%lmaxd, itype, jspin)  = u(1:n_grid_pt, 2, 0:atoms%lmaxd)
            bas1(1:n_grid_pt, 2, 0:atoms%lmaxd, itype, jspin) = du(1:n_grid_pt, 1, 0:atoms%lmaxd)
            bas2(1:n_grid_pt, 2, 0:atoms%lmaxd, itype, jspin) = du(1:n_grid_pt, 2, 0:atoms%lmaxd)

            ! generate radial functions for local orbitals
            IF (atoms%nlo(itype) >= 1) THEN
               CALL radflo(atoms, itype, jspin, enpara%ello0(1, 1, jspin), vr0(:,itype, jspin), &
                           u, du, mpi, usdus, uuilon, duilon, ulouilopn, flo)

               DO ilo = 1, atoms%nlo(itype)
                  bas1(1:n_grid_pt, 2+ilo, atoms%llo(ilo, itype), itype, jspin) = flo(1:n_grid_pt, 1, ilo)
                  bas2(1:n_grid_pt, 2+ilo, atoms%llo(ilo, itype), itype, jspin) = flo(1:n_grid_pt, 2, ilo)
               END DO
            END IF
         END DO
      END DO

      ! the radial functions are normalized
      DO jspin = 1, input%jspins
         DO itype = 1, atoms%ntype
            DO l = 0, atoms%lmax(itype)
               DO i = 1, hybrid%num_radfun_per_l(l, itype)
                  norm = sqrt(intgrf(bas1(:,i, l, itype, jspin)**2 + bas2(:,i, l, itype, jspin)**2, &
                                atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf))
                  bas1(:atoms%jri(itype), i, l, itype, jspin) = bas1(:atoms%jri(itype), i, l, itype, jspin)/norm
                  bas2(:atoms%jri(itype), i, l, itype, jspin) = bas2(:atoms%jri(itype), i, l, itype, jspin)/norm
               END DO
            END DO
         END DO
      END DO
   end subroutine gen_bas_fun

   function calc_selecmat(atoms,hybrid,seleco, selecu) result(selecmat)
      ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
      use m_types
      use m_judft
      implicit NONE

      type(t_atoms),  intent(in) :: atoms
      type(t_hybrid), intent(in) :: hybrid
      LOGICAL, intent(in) :: seleco(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd)
      LOGICAL, intent(in) :: selecu(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd)
      LOGICAL  ::  selecmat(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd, &
                            maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd)
      integer                       :: n1, l1, n2, l2

      ! column-major means left-most index varies the fastest
      do l2=0,atoms%lmaxd
         do n2=1,maxval(hybrid%num_radfun_per_l)
            do l1=0,atoms%lmaxd
               do n1=1,maxval(hybrid%num_radfun_per_l)
                  selecmat(n1,l1,n2,l2) = seleco(n1, l1) .AND. selecu(n2, l2)
               enddo
            enddo
         enddo
      enddo
   end function calc_selecmat

   function calc_radbas_norm(atoms, hybrid, i_basfn, itype, gridf) result(norm)
      use m_types
      implicit NONE
      type(t_atoms), intent(in)  :: atoms
      type(t_hybrid), intent(in) :: hybrid
      integer, intent(in)        :: i_basfn, itype
      real, intent(in)           ::  gridf(:,:)

      real                       :: norm

      norm = SQRT(
                  intgrf(hybrid%radbasfn_mt(:,i_basfn, l, itype)**2, &
                         atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype,&
                         itype, gridf)&
                 )
   end function calc_radbas_norm
END MODULE m_mixedbasis
