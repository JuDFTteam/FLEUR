!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hf_setup

CONTAINS

   SUBROUTINE hf_setup(mpdata, hybinp, input, sym, kpts,  atoms, mpi, noco, cell, oneD, results, jsp, enpara, eig_id_hf, &
                       hybdat, l_real, vr0, eig_irr)
      USE m_types
      USE m_eig66_io
      USE m_util
      USE m_intgrf
      USE m_checkolap
      USE m_hybinp_core
      USE m_gen_wavf
      use m_types_hybdat

      IMPLICIT NONE

      TYPE(t_mpdata), INTENT(inout)   :: mpdata
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_kpts), INTENT(IN)    :: kpts
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
      INTEGER :: nbasfcn, n_dim
      LOGICAL :: l_exist

      ! local arrays

      REAL, ALLOCATABLE :: basprod(:)
      INTEGER              :: degenerat(merge(input%neig*2,input%neig,noco%l_soc) + 1, kpts%nkpt)
      LOGICAL              :: skip_kpt(kpts%nkpt)

      REAL :: zDebug_r(lapw_dim_nbasfcn,input%neig)
      COMPLEX :: zDebug_c(lapw_dim_nbasfcn,input%neig)

      skip_kpt = .FALSE.

      IF (hybinp%l_calhf) THEN
         ! Preparations for HF and hybinp functional calculation
         CALL timestart("gen_bz and gen_wavf")

         allocate(zmat(kpts%nkptf), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation z_c')
         allocate(eig_irr(input%neig, kpts%nkpt), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation eig_irr')
         if(allocated(hybdat%kveclo_eig)) deallocate(hybdat%kveclo_eig)
         allocate(hybdat%kveclo_eig(atoms%nlotot, kpts%nkpt), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%kveclo_eig')
         eig_irr = 0
         hybdat%kveclo_eig = 0

         INQUIRE(file ="z",exist= l_exist)
         IF(l_exist) THEN
            IF (l_real) OPEN(unit=993,file='z',form='unformatted',access='direct',recl=lapw%dim_nbasfcn()*input%neig*8)
            IF (.NOT.l_real) OPEN(unit=993,file='z',form='unformatted',access='direct',recl=lapw%dim_nbasfcn()*input%neig*16)
         END IF

         ! Reading the eig file
         DO nk = 1, kpts%nkpt
            nrec1 = kpts%nkpt*(jsp - 1) + nk
            CALL lapw%init(input, noco, kpts, atoms, sym, nk, cell, sym%zrfs)
            nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
            CALL zMat(nk)%init(l_real, nbasfcn, merge(input%neig*2,input%neig,noco%l_soc))
            CALL read_eig(eig_id_hf, nk, jsp, zmat=zMat(nk))

            IF(l_exist.AND.zmat(1)%l_real) THEN
               READ(993,rec=nk) zDebug_r(:,:)
               zMat(nk)%data_r = 0.0
               zMat(nk)%data_r(:nbasfcn,:input%neig) = zDebug_r(:nbasfcn,:input%neig)
            END IF
            IF(l_exist.AND..NOT.zmat(1)%l_real) THEN
               READ(993,rec=nk) zDebug_c(:,:)
               zMat(nk)%data_c = 0.0
               zMat(nk)%data_c(:nbasfcn,:input%neig) = zDebug_c(:nbasfcn,:input%neig)
            END IF

            eig_irr(:, nk) = results%eig(:, nk, jsp)
            hybinp%ne_eig(nk) = results%neig(nk, jsp)
         END DO

         IF(l_exist) CLOSE(993)

         !Allocate further space
         DO nk = kpts%nkpt + 1, kpts%nkptf
            nbasfcn = zMat(kpts%bkp(nk))%matsize1
            CALL zMat(nk)%init(l_real, nbasfcn, merge(input%neig*2,input%neig,noco%l_soc))
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
         hybinp%nobd(:,jsp) = 0
         DO nk = 1, kpts%nkpt
            DO i = 1, hybinp%ne_eig(nk)
               DO j = i + 1, hybinp%ne_eig(nk)
                  IF (ABS(results%eig(i, nk, jsp) - results%eig(j, nk, jsp)) < 1E-07) THEN !0.015
                     degenerat(i, nk) = degenerat(i, nk) + 1
                  END IF
               END DO
            END DO

            DO i = 1, hybinp%ne_eig(nk)
               IF ((degenerat(i, nk) /= 1) .OR. (degenerat(i, nk) /= 0)) degenerat(i + 1:i + degenerat(i, nk) - 1, nk) = 0
            END DO

            ! set the size of the exchange matrix in the space of the wavefunctions

            hybinp%nbands(nk) = hybinp%bands1
            IF (hybinp%nbands(nk) > hybinp%ne_eig(nk)) THEN
               IF (mpi%irank == 0) THEN
                  WRITE (*, *) ' maximum for hybinp%nbands is', hybinp%ne_eig(nk)
                  WRITE (*, *) ' increase energy window to obtain enough eigenvalues'
                  WRITE (*, *) ' set hybinp%nbands equal to hybinp%ne_eig'
               END IF
               hybinp%nbands(nk) = hybinp%ne_eig(nk)
            END IF

            DO i = hybinp%nbands(nk) - 1, 1, -1
               IF ((degenerat(i, nk) >= 1) .AND. (degenerat(i, nk) + i - 1 /= hybinp%nbands(nk))) THEN
                  hybinp%nbands(nk) = i + degenerat(i, nk) - 1
                  EXIT
               END IF
            END DO

            DO i = 1, hybinp%ne_eig(nk)
               IF (results%w_iks(i, nk, jsp) > 0.0) hybinp%nobd(nk,jsp) = hybinp%nobd(nk,jsp) + 1
            END DO

            IF (hybinp%nobd(nk,jsp) > hybinp%nbands(nk)) THEN
               WRITE (*, *) 'k-point: ', nk
               WRITE (*, *) 'number of bands:          ', hybinp%nbands(nk)
               WRITE (*, *) 'number of occupied bands: ', hybinp%nobd(nk,jsp)
               CALL judft_warn("More occupied bands than total no of bands!?")
               hybinp%nbands(nk) = hybinp%nobd(nk,jsp)
            END IF
            PRINT *, "bands:", nk, hybinp%nobd(nk,jsp), hybinp%nbands(nk), hybinp%ne_eig(nk)
         END DO

         ! spread hybinp%nobd from IBZ to whole BZ
         DO nk = 1, kpts%nkptf
            i = kpts%bkp(nk)
            hybinp%nobd(nk,jsp) = hybinp%nobd(i,jsp)
         END DO

         ! generate eigenvectors z and MT coefficients from the previous iteration at all k-points
         CALL gen_wavf(kpts%nkpt, kpts, sym, atoms, enpara%el0(:, :, jsp), enpara%ello0(:, :, jsp), cell,  &
                       mpdata, hybinp, vr0, hybdat, noco, oneD, mpi, input, jsp, zmat)

         ! generate core wave functions (-> core1/2(jmtd,hybdat%nindxc,0:lmaxc,ntype) )
         CALL corewf(atoms, jsp, input,  vr0, hybdat%lmaxcd, hybdat%maxindxc, mpi, &
                     hybdat%lmaxc, hybdat%nindxc, hybdat%core1, hybdat%core2, hybdat%eig_c)

         ! check olap between core-basis/core-valence/basis-basis
         CALL checkolap(atoms, hybdat, mpdata, hybinp, kpts%nkpt, kpts,  mpi, &
                        input, sym, noco, cell, lapw, jsp)

         ! set up pointer pntgpt

         ! setup dimension of pntgpt
         IF(ALLOCATED(hybdat%pntgptd)) DEALLOCATE(hybdat%pntgptd) ! for spinpolarized systems
         ALLOCATE (hybdat%pntgptd(3))
         hybdat%pntgptd = 0
         DO nk = 1, kpts%nkptf
            CALL lapw%init(input, noco, kpts, atoms, sym, nk, cell, sym%zrfs)
            do n_dim = 1,3
               hybdat%pntgptd(n_dim) = MAXVAL([(ABS(lapw%gvec(n_dim,i,jsp)), i=1, lapw%nv(jsp)), hybdat%pntgptd(n_dim)])
            end do
         END DO

         IF(ALLOCATED(hybdat%pntgpt)) DEALLOCATE(hybdat%pntgpt) ! for spinpolarized systems
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

         allocate(basprod(atoms%jmtd), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation basprod')
         IF(ALLOCATED(hybdat%prodm)) DEALLOCATE(hybdat%prodm)
         allocate(hybdat%prodm(maxval(mpdata%num_radbasfn), hybinp%max_indx_p_1, 0:maxval(hybinp%lcutm1), atoms%ntype), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%prodm')

         call mpdata%init(hybinp, atoms)

         basprod = 0; hybdat%prodm = 0; mpdata%l1 = 0; mpdata%l2 = 0
         mpdata%n1 = 0; mpdata%n2 = 0
         IF(ALLOCATED(hybdat%nindxp1)) DEALLOCATE(hybdat%nindxp1) ! for spinpolarized systems
         ALLOCATE (hybdat%nindxp1(0:maxval(hybinp%lcutm1), atoms%ntype))
         hybdat%nindxp1 = 0
         DO itype = 1, atoms%ntype
            ng = atoms%jri(itype)
            DO l2 = 0, MIN(atoms%lmax(itype), hybinp%lcutwf(itype))
               ll = l2
               DO l1 = 0, ll
                  IF (ABS(l1 - l2) <= hybinp%lcutm1(itype)) THEN
                     DO n2 = 1, mpdata%num_radfun_per_l(l2, itype)
                        nn = mpdata%num_radfun_per_l(l1, itype)
                        IF (l1 == l2) nn = n2
                        DO n1 = 1, nn
                           ! Calculate all basis-function hybdat%products to obtain
                           ! the overlaps with the hybdat%product-basis functions (hybdat%prodm)
                           basprod(:ng) = (hybdat%bas1(:ng, n1, l1, itype)*hybdat%bas1(:ng, n2, l2, itype) + &
                                           hybdat%bas2(:ng, n1, l1, itype)*hybdat%bas2(:ng, n2, l2, itype))/atoms%rmsh(:ng, itype)
                           DO l = ABS(l1 - l2), MIN(hybinp%lcutm1(itype), l1 + l2)
                              IF (MOD(l1 + l2 + l, 2) == 0) THEN
                                 hybdat%nindxp1(l, itype) = hybdat%nindxp1(l, itype) + 1
                                 n = hybdat%nindxp1(l, itype)
                                 mpdata%l1(n,l,itype) = l1
                                 mpdata%l2(n,l,itype) = l2
                                 mpdata%n1(n,l,itype) = n1
                                 mpdata%n2(n,l,itype) = n2
                                 DO i = 1, mpdata%num_radbasfn(l, itype)
                                    hybdat%prodm(i, n, l, itype) = intgrf(basprod(:ng)*mpdata%radbasfn_mt(:ng, i, l, itype), &
                                                                          atoms, itype, hybdat%gridf)
                                 END DO
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
         deallocate(basprod)
         CALL timestop("gen_bz and gen_wavf")

      ELSE IF (hybinp%l_hybrid) THEN ! hybinp%l_calhf is false

         !DO nk = n_start,kpts%nkpt,n_stride
         DO nk = 1, kpts%nkpt, 1
            hybinp%ne_eig(nk) = results%neig(nk, jsp)
            hybinp%nobd(nk,jsp) = COUNT(results%w_iks(:hybinp%ne_eig(nk), nk, jsp) > 0.0)
         END DO

         hybinp%maxlmindx = MAXVAL([(SUM([(mpdata%num_radfun_per_l(l, itype)*(2*l + 1), l=0, atoms%lmax(itype))]), itype=1, atoms%ntype)])
         hybinp%nbands = MIN(hybinp%bands1, input%neig)

      ENDIF ! hybinp%l_calhf

   END SUBROUTINE hf_setup

END MODULE m_hf_setup
