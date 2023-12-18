!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_hf_setup

CONTAINS

   SUBROUTINE hf_setup(mpdata, fi, fmpi,nococonv, results, jsp, enpara, &
                       hybdat, vr0, eig_irr)
      USE m_types
      USE m_constants
      USE m_eig66_io
      USE m_util
      USE m_intgrf
      USE m_checkolap
      USE m_hybrid_core
      USE m_gen_wavf
      use m_types_hybdat

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), INTENT(inout)   :: mpdata
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat

      INTEGER, INTENT(IN)    :: jsp
      REAL, INTENT(IN)    :: vr0(:, :, :)

      REAL, ALLOCATABLE, INTENT(INOUT)   :: eig_irr(:, :)

      ! local type variables
      TYPE(t_lapw)             :: lapw

      ! local scalars
      INTEGER :: ok, nk, nrec1, i, j, l1, l2, ng, itype, n, l, n1, n2, nn
      INTEGER :: nbasfcn, n_dim

      ! local arrays

      REAL, ALLOCATABLE :: basprod(:)
      INTEGER           :: degenerat(merge(fi%input%neig*2,fi%input%neig,fi%noco%l_soc) + 1, fi%kpts%nkpt)


      call timestart("HF_setup")
      call hybdat%set_nobd(fi, results)
      call hybdat%set_nbands(fi, fmpi, results)
      IF (hybdat%l_calhf) THEN
         ! Preparations for HF and hybinp functional calculation
         CALL timestart("gen_bz and gen_wavf")

         ! Reading the eig file
         call timestart("eig stuff")
         DO nk = 1, fi%kpts%nkpt
            nrec1 = fi%kpts%nkpt*(jsp - 1) + nk
            CALL lapw%init(fi%input, fi%noco, nococonv,fi%kpts, fi%atoms, fi%sym, nk, fi%cell)
            nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*fi%atoms%nlotot, lapw%nv(1) + fi%atoms%nlotot, fi%noco%l_noco)

         END DO
         eig_irr = results%eig(:, :, jsp)
         call timestop("eig stuff")


         ! generate eigenvectors z and MT coefficients from the previous iteration at all k-points
         CALL gen_wavf(fi%kpts, fi%sym, fi%atoms, enpara%el0(:, :, jsp), enpara%ello0(:, :, jsp), fi%cell,  &
                       mpdata, vr0, hybdat, fi%noco, nococonv, fmpi, fi%input, jsp)

         ! generate core wave functions (-> core1/2(jmtd,hybdat%nindxc,0:lmaxc,ntype) )
         CALL corewf(fi%atoms, jsp, fi%input,  vr0, hybdat%lmaxcd, hybdat%maxindxc, fmpi, &
                     hybdat%lmaxc, hybdat%nindxc, hybdat%core1, hybdat%core2, hybdat%eig_c)

         ! setup dimension of pntgpt
         IF(ALLOCATED(hybdat%pntgptd)) DEALLOCATE(hybdat%pntgptd) ! for spinpolarized systems
         ALLOCATE (hybdat%pntgptd(3))
         hybdat%pntgptd = 0
         DO nk = 1, fi%kpts%nkptf
            CALL lapw%init(fi%input, fi%noco, nococonv,fi%kpts, fi%atoms, fi%sym, nk, fi%cell)
            do n_dim = 1,3
               hybdat%pntgptd(n_dim) = MAXVAL([(ABS(lapw%gvec(n_dim,i,jsp)), i=1, lapw%nv(jsp)), hybdat%pntgptd(n_dim)])
            end do
         END DO

         IF(ALLOCATED(hybdat%pntgpt)) DEALLOCATE(hybdat%pntgpt) ! for spinpolarized systems
         ALLOCATE (hybdat%pntgpt(-hybdat%pntgptd(1):hybdat%pntgptd(1), -hybdat%pntgptd(2):hybdat%pntgptd(2), &
                                 -hybdat%pntgptd(3):hybdat%pntgptd(3), fi%kpts%nkptf), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation pntgpt')
         hybdat%pntgpt = 0
         DO nk = 1, fi%kpts%nkptf
            CALL lapw%init(fi%input, fi%noco, nococonv,fi%kpts, fi%atoms, fi%sym, nk, fi%cell)
            DO i = 1, lapw%nv(jsp)
               hybdat%pntgpt(lapw%gvec(1,i,jsp), lapw%gvec(2,i,jsp), lapw%gvec(3,i,jsp), nk) = i
            END DO
         END DO

         allocate(basprod(fi%atoms%jmtd), stat=ok, source=0.0)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation basprod')
         IF(ALLOCATED(hybdat%prodm)) DEALLOCATE(hybdat%prodm)
         allocate(hybdat%prodm(maxval(mpdata%num_radbasfn), mpdata%max_indx_p_1, 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), stat=ok)
         IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%prodm')

         hybdat%prodm = 0; mpdata%l1 = 0; mpdata%l2 = 0
         mpdata%n1 = 0; mpdata%n2 = 0
         IF(ALLOCATED(hybdat%nindxp1)) DEALLOCATE(hybdat%nindxp1) ! for spinpolarized systems
         ALLOCATE (hybdat%nindxp1(0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype))
         hybdat%nindxp1 = 0

         !$OMP PARALLEL DO default(none) schedule(dynamic)&
         !$OMP private(itype, ng, l2, l1, n1, l, nn, n, basprod) &
         !$OMP shared(fi, mpdata, hybdat)
         DO itype = 1, fi%atoms%ntype
            ng = fi%atoms%jri(itype)
            DO l2 = 0, MIN(fi%atoms%lmax(itype), fi%hybinp%lcutwf(itype))
               DO l1 = 0, l2
                  IF (ABS(l1 - l2) <= fi%hybinp%lcutm1(itype)) THEN
                     DO n2 = 1, mpdata%num_radfun_per_l(l2, itype)
                        nn = mpdata%num_radfun_per_l(l1, itype)
                        IF (l1 == l2) nn = n2
                        DO n1 = 1, nn
                           ! Calculate all basis-function hybdat%products to obtain
                           ! the overlaps with the hybdat%product-basis functions (hybdat%prodm)
                           basprod(:ng) = (hybdat%bas1(:ng, n1, l1, itype)*hybdat%bas1(:ng, n2, l2, itype) + &
                                           hybdat%bas2(:ng, n1, l1, itype)*hybdat%bas2(:ng, n2, l2, itype))/fi%atoms%rmsh(:ng, itype)
                           DO l = ABS(l1 - l2), MIN(fi%hybinp%lcutm1(itype), l1 + l2)
                              IF (MOD(l1 + l2 + l, 2) == 0) THEN
                                 hybdat%nindxp1(l, itype) = hybdat%nindxp1(l, itype) + 1
                                 n = hybdat%nindxp1(l, itype)
                                 mpdata%l1(n,l,itype) = l1
                                 mpdata%l2(n,l,itype) = l2
                                 mpdata%n1(n,l,itype) = n1
                                 mpdata%n2(n,l,itype) = n2
                                 DO i = 1, mpdata%num_radbasfn(l, itype)
                                    hybdat%prodm(i, n, l, itype) = intgrf(basprod(:ng)*mpdata%radbasfn_mt(:ng, i, l, itype), &
                                                                          fi%atoms, itype, hybdat%gridf)
                                 END DO
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
         deallocate(basprod)
         CALL timestop("gen_bz and gen_wavf")

      ELSE IF (fi%hybinp%l_hybrid) THEN ! hybdat%l_calhf is false
         hybdat%maxlmindx = MAXVAL([(SUM([(mpdata%num_radfun_per_l(l, itype)*(2*l + 1), l=0, fi%atoms%lmax(itype))]), itype=1, fi%atoms%ntype)])
         hybdat%nbands = MIN(fi%hybinp%bands1, fi%input%neig)

      ENDIF ! hybdat%l_calhf

      call timestop("HF_setup")
   END SUBROUTINE hf_setup

END MODULE m_hf_setup
