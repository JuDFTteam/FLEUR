!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hf_setup

CONTAINS

   SUBROUTINE hf_setup(hybrid, input, sym, kpts, DIMENSION, atoms, mpi, noco, cell, oneD, results, jsp, enpara, eig_id_hf, &
                       hybdat, l_real, vr0, eig_irr)
      USE m_types
      USE m_eig66_io
      USE m_util
      USE m_checkolap
      USE m_hybrid_core
      USE m_gen_wavf

      IMPLICIT NONE

      TYPE(t_hybrid), INTENT(INOUT) :: hybrid
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_dimension), INTENT(IN)    :: dimension
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_mpi), INTENT(IN)    :: mpi
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_oneD), INTENT(IN)    :: oneD
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat

      INTEGER, INTENT(IN)    :: jsp, eig_id_hf
      REAL, INTENT(IN)    :: vr0(:, :, :)
      LOGICAL, INTENT(IN)    :: l_real

      REAL, ALLOCATABLE, INTENT(OUT)   :: eig_irr(:, :)

      ! local type variables
      TYPE(t_lapw)             :: lapw
      TYPE(t_mat), ALLOCATABLE :: zmat(:)

      ! local scalars
      INTEGER :: ok, nk, nrec1, i, j, ll, l1, l2, ng, itype, n, l, n1, n2, nn
      INTEGER :: nbasfcn

      ! local arrays

      REAL, ALLOCATABLE :: basprod(:)
      INTEGER              :: degenerat(DIMENSION%neigd2 + 1, kpts%nkpt)
      LOGICAL              :: skip_kpt(kpts%nkpt)

      skip_kpt = .FALSE.

      IF (hybrid%l_calhf) THEN
         ! Preparations for HF and hybrid functional calculation
         CALL timestart("gen_bz and gen_wavf")

         ALLOCATE (zmat(kpts%nkptf), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation z_c')
         ALLOCATE (eig_irr(DIMENSION%neigd2, kpts%nkpt), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation eig_irr')
         ALLOCATE (hybdat%kveclo_eig(atoms%nlotot, kpts%nkpt), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%kveclo_eig')
         eig_irr = 0
         hybdat%kveclo_eig = 0

         ! Reading the eig file
         DO nk = 1, kpts%nkpt
            nrec1 = kpts%nkpt*(jsp - 1) + nk
            CALL lapw%init(input, noco, kpts, atoms, sym, nk, cell, sym%zrfs)
            nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
            CALL zMat(nk)%init(l_real, nbasfcn, dimension%neigd2)
            CALL read_eig(eig_id_hf, nk, jsp, zmat=zMat(nk))
            eig_irr(:, nk) = results%eig(:, nk, jsp)
            hybrid%ne_eig(nk) = results%neig(nk, jsp)
         END DO
         !Allocate further space
         DO nk = kpts%nkpt + 1, kpts%nkptf
            nbasfcn = zMat(kpts%bkp(nk))%matsize1
            CALL zMat(nk)%init(l_real, nbasfcn, dimension%neigd2)
         END DO

         !determine degenerate states at each k-point
         !
         ! degenerat(i) =1  band i  is not degenerat ,
         ! degenerat(i) =j  band i  has j-1 degenart states ( i, i+1, ..., i+j)
         ! degenerat(i) =0  band i  is  degenerat, but is not the lowest band
         !                  of the group of degenerate states
         IF (mpi%irank == 0) THEN
            WRITE (6, *)
            WRITE (6, '(A)') "   k-point      |   number of occupied bands  |   maximal number of bands"
         END IF
         degenerat = 1
         hybrid%nobd = 0
         DO nk = 1, kpts%nkpt
            DO i = 1, hybrid%ne_eig(nk)
               DO j = i + 1, hybrid%ne_eig(nk)
                  IF (ABS(results%eig(i, nk, jsp) - results%eig(j, nk, jsp)) < 1E-07) THEN !0.015
                     degenerat(i, nk) = degenerat(i, nk) + 1
                  END IF
               END DO
            END DO

            DO i = 1, hybrid%ne_eig(nk)
               IF ((degenerat(i, nk) /= 1) .OR. (degenerat(i, nk) /= 0)) degenerat(i + 1:i + degenerat(i, nk) - 1, nk) = 0
            END DO

            ! set the size of the exchange matrix in the space of the wavefunctions

            hybrid%nbands(nk) = hybrid%bands1
            IF (hybrid%nbands(nk) > hybrid%ne_eig(nk)) THEN
               IF (mpi%irank == 0) THEN
                  WRITE (*, *) ' maximum for hybrid%nbands is', hybrid%ne_eig(nk)
                  WRITE (*, *) ' increase energy window to obtain enough eigenvalues'
                  WRITE (*, *) ' set hybrid%nbands equal to hybrid%ne_eig'
               END IF
               hybrid%nbands(nk) = hybrid%ne_eig(nk)
            END IF

            DO i = hybrid%nbands(nk) - 1, 1, -1
               IF ((degenerat(i, nk) >= 1) .AND. (degenerat(i, nk) + i - 1 /= hybrid%nbands(nk))) THEN
                  hybrid%nbands(nk) = i + degenerat(i, nk) - 1
                  EXIT
               END IF
            END DO

            DO i = 1, hybrid%ne_eig(nk)
               IF (results%w_iks(i, nk, jsp) > 0.0) hybrid%nobd(nk) = hybrid%nobd(nk) + 1
            END DO

            IF (hybrid%nobd(nk) > hybrid%nbands(nk)) THEN
               WRITE (*, *) 'k-point: ', nk
               WRITE (*, *) 'number of bands:          ', hybrid%nbands(nk)
               WRITE (*, *) 'number of occupied bands: ', hybrid%nobd(nk)
               CALL judft_warn("More occupied bands than total no of bands!?")
               hybrid%nbands(nk) = hybrid%nobd(nk)
            END IF
            PRINT *, "bands:", nk, hybrid%nobd(nk), hybrid%nbands(nk), hybrid%ne_eig(nk)
         END DO

         ! spread hybrid%nobd from IBZ to whole BZ
         DO nk = 1, kpts%nkptf
            i = kpts%bkp(nk)
            hybrid%nobd(nk) = hybrid%nobd(i)
         END DO

         ! generate eigenvectors z and MT coefficients from the previous iteration at all k-points
         CALL gen_wavf(kpts%nkpt, kpts, sym, atoms, enpara%el0(:, :, jsp), enpara%ello0(:, :, jsp), cell, dimension, &
                       hybrid, vr0, hybdat, noco, oneD, mpi, input, jsp, zmat)

         ! generate core wave functions (-> core1/2(jmtd,hybdat%nindxc,0:lmaxc,ntype) )
         CALL corewf(atoms, jsp, input, DIMENSION, vr0, hybdat%lmaxcd, hybdat%maxindxc, mpi, &
                     hybdat%lmaxc, hybdat%nindxc, hybdat%core1, hybdat%core2, hybdat%eig_c)

         ! check olap between core-basis/core-valence/basis-basis
         CALL checkolap(atoms, hybdat, hybrid, kpts%nkpt, kpts, dimension, mpi, &
                        input, sym, noco, cell, lapw, jsp)

         ! set up pointer pntgpt

         ! setup dimension of pntgpt
         ALLOCATE (hybdat%pntgptd(3))
         hybdat%pntgptd = 0
         DO nk = 1, kpts%nkptf
            CALL lapw%init(input, noco, kpts, atoms, sym, nk, cell, sym%zrfs)
            hybdat%pntgptd(1) = MAXVAL((/(ABS(lapw%k1(i, jsp)), i=1, lapw%nv(jsp)), hybdat%pntgptd(1)/))
            hybdat%pntgptd(2) = MAXVAL((/(ABS(lapw%k2(i, jsp)), i=1, lapw%nv(jsp)), hybdat%pntgptd(2)/))
            hybdat%pntgptd(3) = MAXVAL((/(ABS(lapw%k3(i, jsp)), i=1, lapw%nv(jsp)), hybdat%pntgptd(3)/))
         END DO

         ALLOCATE (hybdat%pntgpt(-hybdat%pntgptd(1):hybdat%pntgptd(1), -hybdat%pntgptd(2):hybdat%pntgptd(2), &
                                 -hybdat%pntgptd(3):hybdat%pntgptd(3), kpts%nkptf), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation pntgpt')
         hybdat%pntgpt = 0
         DO nk = 1, kpts%nkptf
            CALL lapw%init(input, noco, kpts, atoms, sym, nk, cell, sym%zrfs)
            DO i = 1, lapw%nv(jsp)
               hybdat%pntgpt(lapw%gvec(1,i,jsp), lapw%gvec(2,i,jsp), lapw%gvec(3,i,jsp), nk) = i
            END DO
         END DO

         ALLOCATE (basprod(atoms%jmtd), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation basprod')
         ALLOCATE (hybdat%prodm(hybrid%maxindxm1, hybrid%maxindxp1, 0:hybrid%maxlcutm1, atoms%ntype), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%prodm')
         ALLOCATE (hybdat%prod(hybrid%maxindxp1, 0:hybrid%maxlcutm1, atoms%ntype), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%prod')
         basprod = 0; hybdat%prodm = 0; hybdat%prod%l1 = 0; hybdat%prod%l2 = 0
         hybdat%prod%n1 = 0; hybdat%prod%n2 = 0
         ALLOCATE (hybdat%nindxp1(0:hybrid%maxlcutm1, atoms%ntype))
         hybdat%nindxp1 = 0
         DO itype = 1, atoms%ntype
            ng = atoms%jri(itype)
            DO l2 = 0, MIN(atoms%lmax(itype), hybrid%lcutwf(itype))
               ll = l2
               DO l1 = 0, ll
                  IF (ABS(l1 - l2) <= hybrid%lcutm1(itype)) THEN
                     DO n2 = 1, hybrid%nindx(l2, itype)
                        nn = hybrid%nindx(l1, itype)
                        IF (l1 == l2) nn = n2
                        DO n1 = 1, nn
                           ! Calculate all basis-function hybdat%products to obtain
                           ! the overlaps with the hybdat%product-basis functions (hybdat%prodm)
                           basprod(:ng) = (hybdat%bas1(:ng, n1, l1, itype)*hybdat%bas1(:ng, n2, l2, itype) + &
                                           hybdat%bas2(:ng, n1, l1, itype)*hybdat%bas2(:ng, n2, l2, itype))/atoms%rmsh(:ng, itype)
                           DO l = ABS(l1 - l2), MIN(hybrid%lcutm1(itype), l1 + l2)
                              IF (MOD(l1 + l2 + l, 2) == 0) THEN
                                 hybdat%nindxp1(l, itype) = hybdat%nindxp1(l, itype) + 1
                                 n = hybdat%nindxp1(l, itype)
                                 hybdat%prod(n, l, itype)%l1 = l1
                                 hybdat%prod(n, l, itype)%l2 = l2
                                 hybdat%prod(n, l, itype)%n1 = n1
                                 hybdat%prod(n, l, itype)%n2 = n2
                                 DO i = 1, hybrid%nindxm1(l, itype)
                                    hybdat%prodm(i, n, l, itype) = intgrf(basprod(:ng)*hybrid%basm1(:ng, i, l, itype), atoms%jri, &
                                                                          atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)
                                 END DO
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
         DEALLOCATE (basprod)
         CALL timestop("gen_bz and gen_wavf")

      ELSE IF (hybrid%l_hybrid) THEN ! hybrid%l_calhf is false

         !DO nk = n_start,kpts%nkpt,n_stride
         DO nk = 1, kpts%nkpt, 1
            hybrid%ne_eig(nk) = results%neig(nk, jsp)
            hybrid%nobd(nk) = COUNT(results%w_iks(:hybrid%ne_eig(nk), nk, jsp) > 0.0)
         END DO

         hybrid%maxlmindx = MAXVAL((/(SUM((/(hybrid%nindx(l, itype)*(2*l + 1), l=0, atoms%lmax(itype))/)), itype=1, atoms%ntype)/))
         hybrid%nbands = MIN(hybrid%bands1, DIMENSION%neigd)

      ENDIF ! hybrid%l_calhf

   END SUBROUTINE hf_setup

END MODULE m_hf_setup
