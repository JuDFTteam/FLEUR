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
      USE m_util, ONLY: intgrf_init, intgrf, rorderpf
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
      INTEGER                         ::  ikpt, jspin, itype, l1, l2, l, n, igpt, n1, n2, nn, i, j, ng
      INTEGER                         ::  m, nk
      INTEGER                         ::  divconq ! use Divide & Conquer algorithm for array sorting (>0: yes, =0: no)
      REAL                            ::  rdum, rdum1

      ! - local arrays -
      INTEGER                          ::  g(3)
      INTEGER, ALLOCATABLE             ::  ihelp(:)
      INTEGER, ALLOCATABLE             ::  ptr(:)           ! pointer for array sorting
      INTEGER, ALLOCATABLE             ::  unsrt_pgptm(:,:) ! unsorted pointers to g vectors

      REAL                            ::  kvec(3)
      REAL                            ::  bashlp(atoms%jmtd)


      REAL, ALLOCATABLE               ::  olap(:,:), work(:), eig(:), eigv(:,:)
      REAL, ALLOCATABLE               ::  bas1(:,:,:,:,:), bas2(:,:,:,:,:)
      REAL, ALLOCATABLE               ::  basmhlp(:,:,:,:)
      REAL, ALLOCATABLE               ::  gridf(:,:), vr0(:,:,:)
      REAL, ALLOCATABLE               ::  length_kg(:,:) ! length of the vectors k + G

      LOGICAL, ALLOCATABLE            ::  selecmat(:,:,:,:)
      LOGICAL, ALLOCATABLE            ::  seleco(:,:), selecu(:,:)

      CHARACTER, PARAMETER            :: lchar(0:38) = (/'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', &
                                                         'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', &
                                                         'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x'/)


      IF (mpi%irank == 0) WRITE (6, '(//A,I2,A)') '### subroutine: mixedbasis ###'

      IF (xcpot%is_name("exx")) CALL judft_error("EXX is not implemented in this version", calledby='mixedbasis')

      ! Deallocate arrays which might have been allocated in a previous run of this subroutine
      IF (ALLOCATED(mpbasis%ngptm)) deallocate(mpbasis%ngptm)
      IF (ALLOCATED(hybrid%ngptm1)) deallocate(hybrid%ngptm1)
      IF (ALLOCATED(hybrid%nindxm1)) deallocate(hybrid%nindxm1)
      IF (ALLOCATED(hybrid%pgptm)) deallocate(hybrid%pgptm)
      IF (ALLOCATED(hybrid%pgptm1)) deallocate(hybrid%pgptm1)
      IF (ALLOCATED(mpbasis%gptm)) deallocate(mpbasis%gptm)
      IF (ALLOCATED(hybrid%basm1)) deallocate(hybrid%basm1)

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
      call gen_gvec(cell, kpts, mpbasis, hybrid)

      ! construct IR mixed basis set for the representation of the non local exchange elements with cutoff gcutm

      ! first run to determine dimension of pgptm1
      allocate(hybrid%ngptm1(kpts%nkptf))
      hybrid%ngptm1 = 0
      DO igpt = 1, mpbasis%num_gpts()
         g = mpbasis%gptm(:,igpt)
         DO ikpt = 1, kpts%nkptf
            kvec = kpts%bkf(:,ikpt)
            rdum = norm2(MATMUL(kvec + g, cell%bmat))
            IF (rdum <= mpbasis%g_cutoff) THEN
               hybrid%ngptm1(ikpt) = hybrid%ngptm1(ikpt) + 1
            END IF
         END DO
      END DO
      hybrid%maxgptm1 = MAXVAL(hybrid%ngptm1)

      allocate(hybrid%pgptm1(hybrid%maxgptm1, kpts%nkptf))
      hybrid%pgptm1 = 0
      hybrid%ngptm1 = 0

      ! Allocate and initialize arrays needed for G vector ordering
      allocate(unsrt_pgptm(hybrid%maxgptm1, kpts%nkptf))
      allocate(length_kG(hybrid%maxgptm1, kpts%nkptf))
      length_kG = 0
      unsrt_pgptm = 0
      DO igpt = 1, mpbasis%num_gpts()
         g = mpbasis%gptm(:,igpt)
         DO ikpt = 1, kpts%nkptf
            kvec = kpts%bkf(:,ikpt)
            rdum = SUM(MATMUL(kvec + g, cell%bmat)**2)
            IF (rdum <= mpbasis%g_cutoff**2) THEN
               hybrid%ngptm1(ikpt) = hybrid%ngptm1(ikpt) + 1
               unsrt_pgptm(hybrid%ngptm1(ikpt), ikpt) = igpt
               length_kG(hybrid%ngptm1(ikpt), ikpt) = rdum
            END IF
         END DO
      END DO

      ! Sort pointers in array, so that shortest |k+G| comes first
      DO ikpt = 1, kpts%nkptf
         allocate(ptr(hybrid%ngptm1(ikpt)))
         ! Divide and conquer algorithm for arrays > 1000 entries
         divconq = MAX(0, INT(1.443*LOG(0.001*hybrid%ngptm1(ikpt))))
         ! create pointers which correspond to a sorted array
         CALL rorderpf(ptr, length_kG(1:hybrid%ngptm1(ikpt), ikpt), hybrid%ngptm1(ikpt), divconq)
         ! rearrange old pointers
         DO igpt = 1, hybrid%ngptm1(ikpt)
            hybrid%pgptm1(igpt, ikpt) = unsrt_pgptm(ptr(igpt), ikpt)
         END DO
         deallocate(ptr)
      END DO
      deallocate(unsrt_pgptm)
      deallocate(length_kG)

      IF (mpi%irank == 0) THEN
         WRITE (6, '(/A)') 'Mixed basis'
         WRITE (6, '(A,I5)') 'Number of unique G-vectors: ', mpbasis%num_gpts()
         WRITE (6, *)
         WRITE (6, '(3x,A)') 'IR Plane-wave basis with cutoff of gcutm (mpbasis%g_cutoff/2*input%rkmax):'
         WRITE (6, '(5x,A,I5)') 'Maximal number of G-vectors:', maxval(mpbasis%ngptm)
         WRITE (6, *)
         WRITE (6, *)
         WRITE (6, '(3x,A)') 'IR Plane-wave basis for non-local exchange potential:'
         WRITE (6, '(5x,A,I5)') 'Maximal number of G-vectors:', hybrid%maxgptm1
         WRITE (6, *)
      END IF

      ! - - - - - - - - Set up MT product basis for the non-local exchange potential  - - - - - - - - - -

      IF (mpi%irank == 0) THEN
         WRITE (6, '(A)') 'MT product basis for non-local exchange potential:'
         WRITE (6, '(A)') 'Reduction due to overlap (quality of orthonormality, should be < 1.0E-06)'
      END IF

      hybrid%maxlcutm1 = MAXVAL(hybrid%lcutm1)

      allocate(hybrid%nindxm1(0:hybrid%maxlcutm1, atoms%ntype))
      allocate(seleco(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd), selecu(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd))
      allocate(selecmat(maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd, maxval(hybrid%num_radfun_per_l), 0:atoms%lmaxd))
      hybrid%nindxm1 = 0    !!! 01/12/10 jij%M.b.

      ! determine maximal indices of (radial) mixed-basis functions (->nindxm1)
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
            n = 0
            M = 0

            !
            ! valence * valence
            !

            ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
            selecmat = RESHAPE((/((((seleco(n1, l1) .AND. selecu(n2, l2), &
                                     n1=1, maxval(hybrid%num_radfun_per_l)), l1=0, atoms%lmaxd), n2=1, maxval(hybrid%num_radfun_per_l)), l2=0, atoms%lmaxd)/), &
                               [maxval(hybrid%num_radfun_per_l), atoms%lmaxd + 1, maxval(hybrid%num_radfun_per_l), atoms%lmaxd + 1])

            DO l1 = 0, atoms%lmax(itype)
               DO l2 = 0, atoms%lmax(itype)
                  IF (l >= ABS(l1 - l2) .AND. l <= l1 + l2) THEN
                     DO n1 = 1, hybrid%num_radfun_per_l(l1, itype)
                        DO n2 = 1, hybrid%num_radfun_per_l(l2, itype)
                           M = M + 1
                           IF (selecmat(n1, l1, n2, l2)) THEN
                              n = n + 1
                              selecmat(n2, l2, n1, l1) = .FALSE. ! prevent double counting of products (a*b = b*a)
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
            IF (n == 0 .AND. mpi%irank == 0) &
               WRITE (6, '(A)') 'mixedbasis: Warning!  No basis-function product of '//lchar(l)// &
               '-angular momentum defined.'
            hybrid%max_indx_p_1 = MAX(hybrid%max_indx_p_1, M)
            hybrid%nindxm1(l, itype) = n*input%jspins
         END DO
      END DO

      allocate(hybrid%basm1(atoms%jmtd, maxval(hybrid%nindxm1), 0:hybrid%maxlcutm1, atoms%ntype))
      hybrid%basm1 = 0

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
         IF (atoms%ntype > 1 .AND. mpi%irank == 0)&
              &    WRITE (6, '(6X,A,I3)') 'Atom type', itype
         ng = atoms%jri(itype)
         DO l = 0, hybrid%lcutm1(itype)
            n = hybrid%nindxm1(l, itype)
            ! allow for zero product-basis functions for
            ! current l-quantum number
            IF (n == 0) THEN
               IF (mpi%irank == 0) WRITE (6, '(6X,A,'':   0 ->   0'')') lchar(l)
               CYCLE
            END IF

            ! set up the overlap matrix
            allocate(olap(n, n), eigv(n, n), work(3*n), eig(n), ihelp(n))
            ihelp = 1 ! initialize to avoid a segfault
            i = 0

            ! valence*valence

            ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
            selecmat = RESHAPE((/((((seleco(n1, l1) .AND. selecu(n2, l2), &
                                     n1=1, maxval(hybrid%num_radfun_per_l)), l1=0, atoms%lmaxd), n2=1, maxval(hybrid%num_radfun_per_l)), l2=0, atoms%lmaxd)/), &
                               [maxval(hybrid%num_radfun_per_l), atoms%lmaxd + 1, maxval(hybrid%num_radfun_per_l), atoms%lmaxd + 1])

            DO l1 = 0, atoms%lmax(itype)
               DO l2 = 0, atoms%lmax(itype)
                  IF (l < ABS(l1 - l2) .OR. l > l1 + l2) CYCLE

                  DO n1 = 1, hybrid%num_radfun_per_l(l1, itype)
                     DO n2 = 1, hybrid%num_radfun_per_l(l2, itype)

                        IF (selecmat(n1, l1, n2, l2)) THEN
                           DO jspin = 1, input%jspins
                              i = i + 1
                              IF (i > n) call judft_error('got too many product functions', hint='This is a BUG, please report', calledby='mixedbasis')

                              hybrid%basm1(:ng, i, l, itype) &
                                 = (bas1(:ng, n1, l1, itype, jspin) &
                                    *bas1(:ng, n2, l2, itype, jspin) &
                                    + bas2(:ng, n1, l1, itype, jspin) &
                                    *bas2(:ng, n2, l2, itype, jspin))/atoms%rmsh(:ng, itype)

                              !normalize basm1
                              rdum = SQRT(intgrf(hybrid%basm1(:,i, l, itype)**2, &
                                                 atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf))

                              hybrid%basm1(:ng, i, l, itype) = hybrid%basm1(:ng, i, l, itype)/rdum

                           END DO !jspin
                           ! prevent double counting of products (a*b = b*a)
                           selecmat(n2, l2, n1, l1) = .FALSE.
                        END IF

                     END DO !n2
                  END DO !n1

               END DO !l2
            END DO  !l1

            IF (i /= n) call judft_error('counting error for product functions', hint='This is a BUG, please report', calledby='mixedbasis')

            ! In order to get ride of the linear dependencies in the
            ! radial functions basm1 belonging to fixed l and itype
            ! the overlap matrix is diagonalized and those eigenvectors
            ! with a eigenvalue greater then hybrid%tolerance1 are retained

            ! Calculate overlap
            olap = 0
            DO n2 = 1, n
               DO n1 = 1, n2
                  olap(n1, n2) = intgrf(hybrid%basm1(:,n1, l, itype)*hybrid%basm1(:,n2, l, itype), &
                                        atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
                  olap(n2, n1) = olap(n1, n2)
               END DO
            END DO

            ! Diagonalize
            CALL diagonalize(eigv, eig, olap)

            ! Get rid of linear dependencies (eigenvalue <= hybrid%tolerance1)
            nn = 0
            DO i = 1, hybrid%nindxm1(l, itype)
               IF (eig(i) > hybrid%tolerance1) THEN
                  nn = nn + 1
                  ihelp(nn) = i
               END IF
            END DO
            hybrid%nindxm1(l, itype) = nn
            eig = eig(ihelp)
            eigv(:,:) = eigv(:,ihelp)

            DO i = 1, ng
               hybrid%basm1(i, 1:nn, l, itype) = MATMUL(hybrid%basm1(i, 1:n, l, itype), eigv(:,1:nn))/SQRT(eig(:nn))
            END DO

            ! Add constant function to l=0 basis and then do a Gram-Schmidt orthonormalization
            IF (l == 0) THEN

               ! Check if basm1 must be reallocated
               IF (nn + 1 > SIZE(hybrid%basm1, 2)) THEN
                  allocate(basmhlp(atoms%jmtd, nn + 1, 0:hybrid%maxlcutm1, atoms%ntype))
                  basmhlp(:,1:nn, :,:) = hybrid%basm1
                  deallocate(hybrid%basm1)
                  allocate(hybrid%basm1(atoms%jmtd, nn + 1, 0:hybrid%maxlcutm1, atoms%ntype))
                  hybrid%basm1(:,1:nn, :,:) = basmhlp(:,1:nn, :,:)
                  deallocate(basmhlp)
               END IF

               hybrid%basm1(:ng, nn + 1, 0, itype) = atoms%rmsh(:ng, itype)/SQRT(atoms%rmsh(ng, itype)**3/3)
               DO i = nn, 1, -1
                  DO j = i + 1, nn + 1
                     hybrid%basm1(:ng, i, 0, itype) = hybrid%basm1(:ng, i, 0, itype) &
                                                      - intgrf(hybrid%basm1(:ng, i, 0, itype)*hybrid%basm1(:ng, j, 0, itype), &
                                                               atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf) &
                                                      *hybrid%basm1(:ng, j, 0, itype)

                  END DO
                  hybrid%basm1(:ng, i, 0, itype) = hybrid%basm1(:ng, i, 0, itype) &
                                                   /SQRT(intgrf(hybrid%basm1(:ng, i, 0, itype)**2, &
                                                                atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf))
               END DO
               nn = nn + 1
               deallocate(olap)
               allocate(olap(nn, nn))
               hybrid%nindxm1(l, itype) = nn
            END IF

            ! Check orthonormality of product basis
            rdum = 0
            DO i = 1, nn
               DO j = 1, i
                  rdum1 = intgrf(hybrid%basm1(:,i, l, itype)*hybrid%basm1(:,j, l, itype), &
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
               WRITE (6, '(6X,A,I4,'' ->'',I4,''   ('',ES8.1,'' )'')') lchar(l)//':', n, nn, SQRT(rdum)/nn
            END IF

            deallocate(olap, eigv, work, eig, ihelp)

         END DO !l
         IF (mpi%irank == 0) WRITE (6, '(6X,A,I7)') 'Total:', SUM(hybrid%nindxm1(0:hybrid%lcutm1(itype), itype))
      END DO ! itype

      allocate(basmhlp(atoms%jmtd, maxval(hybrid%nindxm1), 0:hybrid%maxlcutm1, atoms%ntype))
      basmhlp(1:atoms%jmtd, 1:maxval(hybrid%nindxm1), 0:hybrid%maxlcutm1, 1:atoms%ntype) &
         = hybrid%basm1(1:atoms%jmtd, 1:maxval(hybrid%nindxm1), 0:hybrid%maxlcutm1, 1:atoms%ntype)
      deallocate(hybrid%basm1)
      allocate(hybrid%basm1(atoms%jmtd, maxval(hybrid%nindxm1), 0:hybrid%maxlcutm1, atoms%ntype))
      hybrid%basm1 = basmhlp

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
         ng = atoms%jri(itype)

         IF (atoms%ntype > 1 .AND. mpi%irank == 0) WRITE (6, '(6X,A,I3)') 'Atom type', itype

         DO l = 0, hybrid%lcutm1(itype)
            ! determine radial function with the largest moment
            ! this function is used to build the linear combinations
            rdum = 0
            DO i = 1, hybrid%nindxm1(l, itype)
               rdum1 = intgrf(atoms%rmsh(:ng, itype)**(l + 1)*hybrid%basm1(:ng, i, l, itype), &
                              atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
               IF (ABS(rdum1) > rdum) THEN
                  n = i
                  rdum = rdum1
               END IF
            END DO

            ! rearrange order of radial functions such that the last function possesses the largest moment
            j = 0
            bashlp(:ng) = hybrid%basm1(:ng, n, l, itype)
            DO i = 1, hybrid%nindxm1(l, itype)
               IF (i == n) CYCLE
               j = j + 1
               hybrid%basm1(:ng, j, l, itype) = hybrid%basm1(:ng, i, l, itype)
            END DO
            hybrid%basm1(:ng, hybrid%nindxm1(l, itype), l, itype) = bashlp(:ng)

         END DO

         DO l = 0, hybrid%lcutm1(itype)
            IF (mpi%irank == 0) WRITE (6, '(6X,A)') lchar(l)//':'

            IF (hybrid%nindxm1(l, itype) == 0) THEN
               IF (mpi%irank == 0) WRITE (6, '(6X,A,'':   0 ->    '')') lchar(l)
               CYCLE
            END IF

            n = hybrid%nindxm1(l, itype)
            DO i = 1, hybrid%nindxm1(l, itype)
               IF (i == n) CYCLE
               ! calculate moment of radial function i
               rdum1 = intgrf(atoms%rmsh(:ng, itype)**(l + 1)*hybrid%basm1(:ng, i, l, itype), &
                              atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

               rdum = intgrf(atoms%rmsh(:ng, itype)**(l + 1)*hybrid%basm1(:ng, n, l, itype), &
                             atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

               bashlp(:ng) = hybrid%basm1(:ng, n, l, itype)

               IF (SQRT(rdum**2 + rdum1**2) <= 1E-06 .AND. mpi%irank == 0) &
                  WRITE (6, *) 'Warning: Norm is smaller thann 1E-06!'

               ! change function n such that n is orthogonal to i
               ! since the functions basm1 have been orthogonal on input
               ! the linear combination does not destroy the orthogonality to the residual functions
               hybrid%basm1(:ng, n, l, itype) = rdum/SQRT(rdum**2 + rdum1**2)*bashlp(:ng) &
                                                + rdum1/SQRT(rdum**2 + rdum1**2)*hybrid%basm1(:ng, i, l, itype)

               ! combine basis function i and n so that they possess no momemt
               hybrid%basm1(:ng, i, l, itype) = rdum1/SQRT(rdum**2 + rdum1**2)*bashlp(:ng) &
                                                - rdum/SQRT(rdum**2 + rdum1**2)*hybrid%basm1(:ng, i, l, itype)

               rdum1 = intgrf(atoms%rmsh(:ng, itype)**(l + 1)*hybrid%basm1(:ng, i, l, itype), &
                              atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

               IF (rdum1 > 1E-10) call judft_error('moment of radial function does not vanish', calledby='mixedbasis')

               IF (mpi%irank == 0) WRITE (6, '(6x,I4,'' ->  '',ES8.1)') i, rdum1
            END DO

            ! test orthogonality
            rdum = 0
            DO i = 1, hybrid%nindxm1(l, itype)
               DO j = 1, hybrid%nindxm1(l, itype)
                  rdum1 = intgrf(hybrid%basm1(:ng, i, l, itype)*hybrid%basm1(:ng, j, l, itype), &
                                 atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
                  IF (i /= j) THEN
                     rdum = rdum + rdum1
                  ELSE
                     rdum = rdum + ABS(1 - rdum1)
                  END IF
               END DO
            END DO
            IF (mpi%irank == 0) &
               WRITE (6, '(6x,I4,'' ->'',f10.5,''   ('',ES8.1,'' )'')') n, &
               intgrf(atoms%rmsh(:ng, itype)**(l + 1)*hybrid%basm1(:ng, n, l, itype), atoms%jri, &
                      atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf), rdum
         END DO
      END DO

      DO itype = 1, atoms%ntype
         IF (ANY(hybrid%nindxm1(0:hybrid%lcutm1(itype), itype) == 0)) call judft_error('any hybrid%nindxm1 eq 0', calledby='mixedbasis')
      END DO

      !count basis functions
      hybrid%nbasp = 0
      DO itype = 1, atoms%ntype
         DO i = 1, atoms%neq(itype)
            DO l = 0, hybrid%lcutm1(itype)
               DO M = -l, l
                  DO j = 1, hybrid%nindxm1(l, itype)
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
      USE m_util,   ONLY: intgrf
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
      INTEGER :: ng, noded, nodem
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
         ng = atoms%jri(itype) ! number of radial gridpoints
         DO jspin = 1, input%jspins
            DO l = 0, atoms%lmax(itype)
               CALL radfun(l, itype, jspin, enpara%el0(l, itype, jspin), vr0(:,itype, jspin), atoms, &
                           u(:,:,l), du(:,:,l), usdus, nodem, noded, wronk)
            END DO
            bas1(1:ng, 1, 0:atoms%lmaxd, itype, jspin)  = u(1:ng, 1, 0:atoms%lmaxd)
            bas2(1:ng, 1, 0:atoms%lmaxd, itype, jspin)  = u(1:ng, 2, 0:atoms%lmaxd)
            bas1(1:ng, 2, 0:atoms%lmaxd, itype, jspin) = du(1:ng, 1, 0:atoms%lmaxd)
            bas2(1:ng, 2, 0:atoms%lmaxd, itype, jspin) = du(1:ng, 2, 0:atoms%lmaxd)

            ! generate radial functions for local orbitals
            IF (atoms%nlo(itype) >= 1) THEN
               CALL radflo(atoms, itype, jspin, enpara%ello0(1, 1, jspin), vr0(:,itype, jspin), &
                           u, du, mpi, usdus, uuilon, duilon, ulouilopn, flo)

               DO ilo = 1, atoms%nlo(itype)
                  bas1(1:ng, 2+ilo, atoms%llo(ilo, itype), itype, jspin) = flo(1:ng, 1, ilo)
                  bas2(1:ng, 2+ilo, atoms%llo(ilo, itype), itype, jspin) = flo(1:ng, 2, ilo)
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

   subroutine gen_gvec(cell, kpts, mpbasis, hybrid)
      use m_types
      USE m_util, ONLY: intgrf_init, intgrf, rorderpf
      implicit NONE
      type(t_cell), intent(in)       :: cell
      type(t_kpts), intent(in)       :: kpts
      TYPE(t_mpbasis), intent(inout) :: mpbasis
      type(t_hybrid), intent(inout)  :: hybrid


      integer :: i, n, n1, n2, divconq
      integer :: x, y, z, ikpt, igpt
      integer :: g(3)
      real    :: longest_k, kvec(3)
      logical :: l_found_new_gpt, l_found_kg_in_sphere

      INTEGER, ALLOCATABLE            ::  unsrt_pgptm(:,:) ! unsorted pointers to g vectors
      REAL, ALLOCATABLE               ::  length_kg(:,:) ! length of the vectors k + G
      INTEGER, ALLOCATABLE            ::  ptr(:)

      allocate(mpbasis%ngptm(kpts%nkptf))

      mpbasis%ngptm = 0
      i = 0
      n = -1

      longest_k = MAXVAL([(norm2(MATMUL(kpts%bkf(:,ikpt), cell%bmat)), ikpt=1, kpts%nkptf)])

      ! a first run for the determination of the dimensions of the fields gptm,pgptm

      DO
         n = n + 1
         l_found_new_gpt = .FALSE.
         DO x = -n, n
            n1 = n - ABS(x)
            DO y = -n1, n1
               n2 = n1 - ABS(y)
               DO z = -n2, n2, MAX(2*n2, 1)
                  g = [x, y, z]
                  IF ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpbasis%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  DO ikpt = 1, kpts%nkptf
                     IF (norm2(MATMUL(kpts%bkf(:,ikpt) + g, cell%bmat)) <= mpbasis%g_cutoff) THEN
                        IF (.NOT. l_found_kg_in_sphere) THEN
                           i = i + 1
                           l_found_kg_in_sphere = .TRUE.
                        END IF

                        mpbasis%ngptm(ikpt) = mpbasis%ngptm(ikpt) + 1
                        l_found_new_gpt = .TRUE.
                     END IF
                  END DO ! k-loop
               END DO
            END DO
         END DO
         IF (.NOT. l_found_new_gpt) EXIT
      END DO

      allocate(mpbasis%gptm(3,i)) ! i = gptmd
      allocate(hybrid%pgptm(maxval(mpbasis%ngptm), kpts%nkptf))

      ! Allocate and initialize arrays needed for G vector ordering
      allocate(unsrt_pgptm(maxval(mpbasis%ngptm), kpts%nkptf))
      allocate(length_kG(maxval(mpbasis%ngptm), kpts%nkptf))

      mpbasis%gptm = 0
      hybrid%pgptm = 0
      mpbasis%ngptm = 0

      i = 0
      n = -1

      length_kG = 0
      unsrt_pgptm = 0

      DO
         n = n + 1
         l_found_new_gpt = .FALSE.
         DO x = -n, n
            n1 = n - ABS(x)
            DO y = -n1, n1
               n2 = n1 - ABS(y)
               DO z = -n2, n2, MAX(2*n2, 1)
                  g = [x, y, z]
                  IF ((norm2(MATMUL(g, cell%bmat)) - longest_k) > mpbasis%g_cutoff) CYCLE
                  l_found_kg_in_sphere = .FALSE.
                  DO ikpt = 1, kpts%nkptf
                     kvec = kpts%bkf(:,ikpt)

                     IF (norm2(MATMUL(kvec + g, cell%bmat)) <= mpbasis%g_cutoff) THEN
                        IF (.NOT. l_found_kg_in_sphere) THEN
                           i = i + 1
                           mpbasis%gptm(:,i) = g
                           l_found_kg_in_sphere = .TRUE.
                        END IF

                        mpbasis%ngptm(ikpt) = mpbasis%ngptm(ikpt) + 1
                        l_found_new_gpt = .TRUE.

                        ! Position of the vector is saved as pointer
                        unsrt_pgptm(mpbasis%ngptm(ikpt), ikpt) = i
                        ! Save length of vector k + G for array sorting
                        length_kG(mpbasis%ngptm(ikpt), ikpt) = norm2(MATMUL(kvec + g, cell%bmat))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         IF (.NOT. l_found_new_gpt) EXIT
      END DO

      ! Sort pointers in array, so that shortest |k+G| comes first
      DO ikpt = 1, kpts%nkptf
         allocate(ptr(mpbasis%ngptm(ikpt)))
         ! Divide and conquer algorithm for arrays > 1000 entries
         divconq = MAX(0, INT(1.443*LOG(0.001*mpbasis%ngptm(ikpt))))
         ! create pointers which correspond to a sorted array
         CALL rorderpf(ptr, length_kG(1:mpbasis%ngptm(ikpt), ikpt), mpbasis%ngptm(ikpt), divconq)
         ! rearrange old pointers
         DO igpt = 1, mpbasis%ngptm(ikpt)
            hybrid%pgptm(igpt, ikpt) = unsrt_pgptm(ptr(igpt), ikpt)
         END DO
         deallocate(ptr)
      END DO
   end subroutine gen_gvec
END MODULE m_mixedbasis
