!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_nonsph
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC hsmt_nonsph

CONTAINS
   SUBROUTINE hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,igSpinPr,igSpin,chi,noco,nococonv,cell,lapw,td,fjgj,hmat,set0,lapwq,fjgjq)
      USE m_hsmt_fjgj
      USE m_types
      USE m_hsmt_ab
#ifdef _OPENACC
      USE cublas
#define CPP_zgemm cublaszgemm
#define CPP_zherk cublaszherk
#define CPP_data_c data_c
#else
#define CPP_zgemm zgemm
#define CPP_zherk zherk
#define CPP_data_c hmat%data_c
#endif

      TYPE(t_mpi),      INTENT(IN) :: fmpi
      TYPE(t_sym),      INTENT(IN) :: sym
      TYPE(t_noco),     INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_cell),     INTENT(IN) :: cell
      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_lapw),     INTENT(IN) :: lapw
      TYPE(t_tlmplm),   INTENT(IN) :: td
      TYPE(t_fjgj),     INTENT(IN) :: fjgj

      INTEGER,          INTENT(IN) :: n, ilSpinPr, ilSpin, igSpinPr, igSpin
      COMPLEX,          INTENT(IN) :: chi
      LOGICAL,          INTENT(IN) :: set0  !if true, initialize the hmat matrix with zeros

      CLASS(t_mat),INTENT(INOUT)     ::hmat

      TYPE(t_lapw), OPTIONAL, INTENT(IN) :: lapwq ! Additional set of lapw, in case
      TYPE(t_fjgj), OPTIONAL, INTENT(IN) :: fjgjq ! the left and right ones differ.

      INTEGER :: nn, na, ab_size, l, ll, size_ab_select
      INTEGER :: size_data_c, size_ab, size_ab2 !these data-dimensions are not available so easily in openacc, hence we store them
      INTEGER :: ikGPr, ikG
      REAL    :: rchi
      COMPLEX :: cchi
      LOGICAL :: l_samelapw

      COMPLEX, ALLOCATABLE :: ab1(:,:),ab_select(:,:)
      COMPLEX, ALLOCATABLE :: abCoeffs(:,:), ab2(:,:), h_loc(:,:), data_c(:,:)

      COMPLEX, ALLOCATABLE :: abCoeffsPr(:,:)
      TYPE(t_lapw) :: lapwPr
      TYPE(t_fjgj) :: fjgjPr

      CALL timestart("non-spherical setup")

      l_samelapw = .FALSE.
      IF (.NOT.PRESENT(lapwq)) l_samelapw = .TRUE.
      IF (.NOT.l_samelapw) THEN
         lapwPr = lapwq
         fjgjPr = fjgjq
      ELSE
         lapwPr = lapw
         fjgjPr = fjgj
      END IF

      size_ab = maxval(lapw%nv)

      IF (fmpi%n_size==1) Then
         size_ab_select=size_ab
      ELSE
         size_ab_select=lapwPr%num_local_cols(igSpinPr)
      END IF

      ALLOCATE(ab_select(size_ab_select, 2 * atoms%lmaxd * (atoms%lmaxd + 2) + 2))
      ALLOCATE(abCoeffs(2 * atoms%lmaxd * (atoms%lmaxd + 2) + 2, MAXVAL(lapw%nv)),&
             & ab1(size_ab, 2 * atoms%lmaxd * (atoms%lmaxd + 2) + 2))
      ! TODO: Check, whether this is necessary or shifting to
      !       max(MAXVAL(lapwq%nv),MAXVAL(lapw%nv)) in abCoeffs is also enough.
      ALLOCATE(abCoeffsPr(2 * atoms%lmaxd * (atoms%lmaxd + 2) + 2, MAXVAL(lapwPr%nv)))

      IF (igSpinPr.NE.igSpin) THEN
         ALLOCATE(ab2(lapwPr%nv(igSpinPr), 2 * atoms%lmaxd * (atoms%lmaxd + 2) + 2))
         size_ab2 = lapwPr%nv(igSpinPr)
      ELSE
         ALLOCATE(ab2(1,1))
         size_ab2 = 1
      END IF

#ifndef _OPENACC
      IF (hmat%l_real) THEN
         IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
            DEALLOCATE(hmat%data_c)
            ALLOCATE(hmat%data_c(SIZE(hmat%data_r, 1), SIZE(hmat%data_r, 2)))
         END IF
         !$OMP PARALLEL DO DEFAULT(shared)
         DO l = 1, size(hmat%data_c, 2)
            hmat%data_c(:,l) = 0.0
         END DO
         !$OMP END PARALLEL DO
      END IF
      size_data_c = size(hmat%data_c, 1)
#else
      IF (hmat%l_real) THEN
         ALLOCATE(data_c(SIZE(hmat%data_r, 1), SIZE(hmat%data_r, 2)))
         size_data_c = size(data_c, 1)
      ELSE
         ALLOCATE(data_c(SIZE(hmat%data_c, 1), SIZE(hmat%data_c, 2)))
         size_data_c = size(data_c, 1)
      END IF
#endif

      ALLOCATE(h_loc(SIZE(td%h_loc_nonsph, 1), SIZE(td%h_loc_nonsph, 1)))
      h_loc = td%h_loc_nonsph(0:, 0:, n, ilSpinPr, ilSpin)

#ifdef _OPENACC
      !$acc enter data create(ab2,ab1,abCoeffs,abCoeffsPr,data_c,ab_select)copyin(h_loc)
      !$acc kernels present(data_c) default(none)
      data_c(:, :)=0.0
      !$acc end kernels
#endif

      DO nn = 1,atoms%neq(n)
         na = SUM(atoms%neq(:n-1)) + nn
         IF ((sym%invsat(na)==0) .OR. (sym%invsat(na)==1)) THEN
            rchi = MERGE(REAL(chi), REAL(chi)*2, (sym%invsat(na)==0))
            cchi = MERGE(chi, chi*2, (sym%invsat(na)==0))

            ! abCoeffs for \sigma_{\alpha} and \sigma_{g}
            ! Denoted in comments as a
            ! [local spin primed -> '; global spin primed -> pr]
            CALL timestart("hsmt_ab_1")
            CALL hsmt_ab(sym, atoms, noco, nococonv, ilSpin, igSpin, n, na, cell, &
                       & lapw, fjgj, abCoeffs, ab_size, .TRUE.)
            CALL timestop("hsmt_ab_1")

            IF (l_samelapw.AND.(ilSpinPr==ilSpin)) THEN
               !!$acc update device(ab)
               !$acc host_data use_device(abCoeffs,ab1,h_loc)
               CALL CPP_zgemm("C", "N", lapw%nv(igSpin), ab_size, ab_size, cmplx(1.0, 0.0), &
                            & abCoeffs, SIZE(abCoeffs, 1), h_loc, size(td%h_loc_nonsph, 1), &
                            & cmplx(0.0, 0.0), ab1, size_ab)
               !$acc end host_data
            ELSE ! Needed, because t^H .NE. t!
               !!$acc update device(ab)
               !$acc host_data use_device(abCoeffs,ab1,h_loc)
               CALL CPP_zgemm("C", "C", lapw%nv(igSpin), ab_size, ab_size, cmplx(1.0, 0.0), &
                            & abCoeffs, SIZE(abCoeffs, 1), h_loc, size(td%h_loc_nonsph, 1), &
                            & cmplx(0.0, 0.0), ab1, size_ab)
               !$acc end host_data
            END IF

            ! ab1 = MATMUL(TRANSPOSE(abCoeffs(:ab_size,:lapw%nv(igSpin))),h_loc(:ab_size,:ab_size,n,ilSpin))
            ! In locally diagonal case:
            ! ab1 = a^H * L (lower triangular matrix from Cholesky decomposition)
            ! Locally offdiagonal case:
            ! ab1 = a^H * t (potential matrix in lmp lm etc.)
            ! .NOT.l_samelapw:
            ! ab1 = a^H * t^H

            ! Of these ab1 coeffs only a part is needed in case of MPI parallelism
            !$acc kernels default(none) present(ab_select,ab1)copyin(fmpi)
            IF (fmpi%n_size>1) THEN
               ab_select(:, :) = ab1(fmpi%n_rank+1:lapw%nv(igSpin):fmpi%n_size, :)
            ELSE
               ab_select(:, :) = ab1(:, :) !All of ab1 needed
            END IF
            !$acc end kernels

            IF (igSpinPr==igSpin) THEN
               IF (ilSpinPr==ilSpin) THEN
                  IF (l_samelapw) THEN
                     IF (fmpi%n_size==1) THEN !use z-herk trick on single PE
                        !$acc host_data use_device(data_c,ab1)
                        IF (set0 .and. nn == 1) THEN
                           !CPP_data_c = CMPLX(0.0,0.0)
                           CALL CPP_zherk("U", "N", lapw%nv(igSpinPr), ab_size, Rchi, &
                                        & ab1, size_ab, 0.0, CPP_data_c, size_data_c)
                        ELSE
                           CALL CPP_zherk("U", "N", lapw%nv(igSpinPr), ab_size, Rchi, &
                                        & ab1, size_ab, 1.0, CPP_data_c, size_data_c)
                        END IF
                        !$acc end host_data
                        ! conjgsolve:
                        ! data_c += Rchi * a^H * H * a
                        ! [only upper triangle]
                     ELSE ! zgemm case
                        !$acc host_data use_device(data_c,ab1,ab_select)
                        IF (set0 .and. nn == 1) THEN
                           !CPP_data_c = CMPLX(0.0,0.0)
                           CALL CPP_zgemm("N", "C", lapw%nv(igSpinPr), size_ab_select, ab_size, cchi, &
                                        & ab1, size_ab, ab_select, lapw%num_local_cols(igSpinPr), &
                                        & CMPLX(0.0, 0.0), CPP_data_c, size_data_c)
                        ELSE
                           CALL CPP_zgemm("N", "C", lapw%nv(igSpinPr), size_ab_select, ab_size, cchi, &
                                        & ab1, size_ab, ab_select, lapw%num_local_cols(igSpinPr), &
                                        & CMPLX(1.0, 0.0), CPP_data_c, size_data_c)
                        END IF
                        !$acc end host_data
                        ! conjgsolve:
                        ! ab_select = a^H * L
                        ! ab1 = a^H * L
                        ! data_c += cchi * ab1 * abselect^H
                        !         = cchi * a^H * H * a
                     END IF
                  ELSE ! Case for additional q on left vector.
                     CALL timestart("hsmt_ab_2")
                     CALL hsmt_ab(sym, atoms, noco, nococonv, ilSpin, igSpin, n, na, cell, &
                                & lapwPr, fjgjPr, abCoeffsPr, ab_size, .TRUE.)
                     !!$acc update device (abCoeffsPr)
                     CALL timestop("hsmt_ab_2")

                     !$acc host_data use_device(abCoeffsPr,data_c,ab1,ab_select)
                     IF (set0 .and. nn == 1) THEN
                        CALL CPP_zgemm("C", "C", lapwPr%nv(igSpin), size_ab_select, ab_size, chi, &
                                     & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                     & CMPLX(0.0, 0.0), CPP_data_c, SIZE_data_c)
                     ELSE
                        CALL CPP_zgemm("C", "C", lapwPr%nv(igSpin), size_ab_select, ab_size, chi, &
                                     & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                     & CMPLX(1.0, 0.0), CPP_data_c, SIZE_data_c)
                     END IF
                     !$acc end host_data
                     ! data_c += chi * aq * abselect^H
                     !         = chi * aq^H * t * a
                  END IF
               ELSE !This is the case of a local off-diagonal contribution.
                  !It is not Hermitian, so we NEED to use zgemm CALL

                  ! abCoeffs for \sigma_{\alpha}^{'} and \sigma_{g}
                  CALL timestart("hsmt_ab_3")
                  CALL hsmt_ab(sym, atoms, noco, nococonv, ilSpinPr, igSpin, n, na, cell, &
                             & lapwPr, fjgjPr, abCoeffsPr, ab_size, .TRUE.)
                  !!$acc update device(abCoeffsPr)
                  CALL timestop("hsmt_ab_3")

                  !$acc host_data use_device(abCoeffs,data_c,ab1,ab_select)
                  IF (set0 .and. nn == 1) THEN
                     !CPP_data_c = CMPLX(0.0,0.0)
                     CALL CPP_zgemm("C", "C", lapwPr%nv(igSpinPr), size_ab_select, ab_size, chi, &
                                  & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                  & CMPLX(0.0, 0.0), CPP_data_c, SIZE_data_c)
                  ELSE
                     CALL CPP_zgemm("C", "C", lapwPr%nv(igSpinPr), size_ab_select, ab_size, chi, &
                                  & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                  & CMPLX(1.0, 0.0), CPP_data_c, SIZE_data_c)
                  END IF
                  !$acc end host_data
                  ! conjgsolve:
                  ! ab_select = a^H * t
                  ! abCoeffs = a'
                  ! data_c += chi * abCoeffs^H * ab_select^H
                  !         = chi * a'^H * t * a
                  ! .NOT.l_samelapw:
                  ! ab_select = a^H * t^H
                  ! abCoeffs = aq'
                  ! data_c += chi * abCoeffs^H * ab_select^H
                  !         = chi * aq'^H * t * a
               END IF
            ELSE  !here the l_ss off-diagonal part starts
               !Second set of abCoeffs is needed
               ! abCoeffs for \sigma_{\alpha}^{'} and \sigma_{g}^{'}
               CALL timestart("hsmt_ab_4")
               CALL hsmt_ab(sym, atoms, noco, nococonv, ilSpinPr, igSpinPr, n, na, cell, &
                          & lapwPr, fjgjPr, abCoeffsPr, ab_size, .TRUE.)
               CALL timestop("hsmt_ab_4")
               IF (ilSpinPr==ilSpin) THEN
                  IF (l_samelapw) THEN
                     !!$acc update device (abCoeffs)
                     !$acc host_data use_device(abCoeffs,h_loc,ab2)
                     CALL CPP_zgemm("C", "N", lapwPr%nv(igSpinPr), ab_size, ab_size, CMPLX(1.0, 0.0), &
                                  & abCoeffsPr, SIZE(abCoeffsPr, 1), h_loc, size(td%h_loc_nonsph, 1), &
                                  & CMPLX(0.0, 0.0), ab2, size_ab2)
                     !$acc end host_data
                     !Multiply for Hamiltonian
                     !$acc host_data use_device(ab2,ab1,data_c,ab_select)
                     CALL CPP_zgemm("N", "C", lapwPr%nv(igSpinPr), lapwPr%num_local_cols(igSpin), ab_size, chi, &
                                  & ab2, size_ab2, ab_select, size_ab_select, &
                                  & CMPLX(1.0, 0.0), CPP_data_c, size_data_c)
                     !$acc end host_data
                     ! conjgsolve:
                     ! ab2 = aPr'^H * L
                     ! ab_select = a^H * L
                     ! data_c += chi * ab2 * ab_select^H
                     !         = chi * aPr'^H * H * a
                  ELSE
                     !$acc host_data use_device(abCoeffs,data_c,ab_select)
                     IF (set0 .AND. nn == 1) THEN
                        !CPP_data_c = CMPLX(0.0,0.0)
                        CALL CPP_zgemm("C", "C", lapwPr%nv(igSpinPr), lapwPr%num_local_cols(igSpin), ab_size, cchi, &
                                     & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                     & CMPLX(0.0, 0.0), CPP_data_c, SIZE_data_c)
                     ELSE
                        CALL CPP_zgemm("C", "C", lapwPr%nv(igSpinPr), lapwPr%num_local_cols(igSpin), ab_size, cchi, &
                                     & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                     CMPLX(1.0, 0.0), CPP_data_c, SIZE_data_c)
                     END IF
                     !$acc end host_data
                     ! abCoeffs = aqPr'
                     ! ab_select = a^H * t^H
                     ! data_c += cchi * abCoeffs^H *  abselect^H
                     !         = cchi * aqPr'^H * t * a
                  END IF
               ELSE
                  !$acc host_data use_device(abCoeffs,ab1,data_c,ab_select)
                  IF (set0 .AND. nn == 1) THEN
                     !CPP_data_c = CMPLX(0.0,0.0)
                     CALL CPP_zgemm("C", "C", lapwPr%nv(igSpinPr), lapwPr%num_local_cols(igSpin), ab_size, cchi, &
                                  & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                  & CMPLX(0.0, 0.0), CPP_data_c, SIZE_data_c)
                  ELSE
                     CALL CPP_zgemm("C", "C", lapwPr%nv(igSpinPr), lapwPr%num_local_cols(igSpin), ab_size, cchi, &
                                  & abCoeffsPr, SIZE(abCoeffsPr, 1), ab_select, size_ab_select, &
                                  CMPLX(1.0, 0.0), CPP_data_c, SIZE_data_c)
                  END IF
                  !$acc end host_data
                  ! conjgsolve:
                  ! ab_select = a^H * t
                  ! abCoeffs = aPr'
                  ! data_c += chi * abCoeffs^H * ab_select^H
                  !         = chi * aPr'^H * t * a
                  ! .NOT.l_samelapw:
                  ! ab_select = a^H * t^H
                  ! abCoeffs = aqPr'
                  ! data_c += chi * abCoeffs^H * ab_select^H
                  !         = chi * aqPr'^H * t * a
               END IF
            END IF
         END IF
      END DO

#ifdef _OPENACC
      IF (hmat%l_real) THEN
         !$acc kernels present(hmat,hmat%data_r,data_c) default(none)
         hmat%data_r = hmat%data_r + real(data_c)
         !$acc end kernels
      ELSE
         !$acc kernels present(hmat,hmat%data_c,data_c) default(none)
         hmat%data_c = hmat%data_c + data_c
         !$acc end kernels
      END IF
#else
      IF (hmat%l_real) THEN
         !$OMP PARALLEL DO DEFAULT(shared)
         DO l = 1, size(hmat%data_c, 2)
            hmat%data_r(:, l) = hmat%data_r(:, l) + REAL(hmat%data_c(:, l))
         END DO
         !$OMP END PARALLEL DO
      END IF
#endif

      !$acc exit data delete(ab2,ab1,abCoeffs,abCoeffsPr,data_c,ab_select,h_loc)

      CALL timestop("non-spherical setup")
   END SUBROUTINE hsmt_nonsph

END MODULE m_hsmt_nonsph
