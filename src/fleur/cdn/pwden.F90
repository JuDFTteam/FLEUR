!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_pwden
   USE m_types_dos
   USE m_constants
   USE m_forceb8
   USE m_pwint
   USE m_juDFT
   USE m_types_fftGrid
   USE m_types_fft
   USE m_selecFFT
   USE m_fft_interface
#ifdef _OPENACC
   USE openacc
#endif
   !$ USE omp_lib
   USE m_types_atoms
   USE m_types_banddos
   USE m_types_cell
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_mat
   USE m_types_mpi
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_potden
   USE m_types_misc
   USE m_types_stars
   USE m_types_sym
   implicit none
   PRIVATE
   PUBLIC :: pwden
CONTAINS
   SUBROUTINE pwden(stars, kpts, banddos,   input, fmpi, noco, nococonv, cell, atoms, sym, &
                    ikpt, jspin, lapw, ne, ev_list, we, eig, den, results, f_b8, zMat, dos, q_dfpt, lapwq, we1, zMat1, iDir, &
                    lapwmq,zMat1m)
      !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !     In this subroutine the star function expansion coefficients of
      !     the plane wave charge density is determined.
      !
      !     This subroutine is called for each k-point and each spin.
      !
      !
      !     INPUT:    eigen vectors
      !               reciprocal lattice information
      !               Brillouine zone sampling
      !               FFT information
      !
      !     OUTPUT:   den%pw(s)
      !
      !                                      Stefan Bl"ugel, JRCAT, Feb. 1997
      !                                      Gustav Bihlmayer, UniWien
      !
      !     In non-collinear calculations the density becomes a hermitian 2x2
      !     matrix. This subroutine generates this density matrix in the
      !     interstitial region. The diagonal elements of this matrix
      !     (n_11 & n_22) are stored in den%pw, while the real and imaginary part
      !     of the off-diagonal element are store in den%pw(:,3).
      !
      !     Philipp Kurz 99/07
      !
      !     Subtroutine was refactored in 2020. GM.
      !
      !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      IMPLICIT NONE

      TYPE(t_lapw), INTENT(IN)       :: lapw
      TYPE(t_mpi), INTENT(IN)        :: fmpi

      TYPE(t_banddos), INTENT(IN)    :: banddos
      TYPE(t_input), INTENT(IN)      :: input
      TYPE(t_noco), INTENT(IN)       :: noco
      TYPE(t_nococonv),INTENT(IN)    :: nococonv
      TYPE(t_sym), INTENT(IN)        :: sym
      TYPE(t_stars), INTENT(IN)      :: stars
      TYPE(t_cell), INTENT(IN)       :: cell
      TYPE(t_kpts), INTENT(IN)       :: kpts
      TYPE(t_atoms), INTENT(IN)      :: atoms
      TYPE(t_mat), INTENT(IN)        :: zMat
      TYPE(t_potden), INTENT(INOUT)  :: den
      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_dos), INTENT(INOUT)     :: dos

      REAL, INTENT(IN)       :: we(:)  !(nobd)
      REAL, INTENT(IN)       :: eig(:) !(input%neig)
      INTEGER, INTENT(IN)    :: ne
      INTEGER, INTENT(IN)    :: ev_list(ne)
      INTEGER, INTENT(IN)    :: ikpt
      INTEGER, INTENT(IN)    :: jspin
      COMPLEX, INTENT(INOUT) ::  f_b8(3, atoms%ntype)

      REAL,         OPTIONAL, INTENT(IN) :: q_dfpt(3), we1(:)
      TYPE(t_mat),  OPTIONAL, INTENT(IN) :: zMat1
      TYPE(t_lapw), OPTIONAL, INTENT(IN) :: lapwq
      INTEGER,      OPTIONAL, INTENT(IN) :: iDir

      TYPE(t_lapw), OPTIONAL, INTENT(IN) :: lapwmq
      TYPE(t_mat),  OPTIONAL, INTENT(IN) :: zMat1m

      ! local variables
      TYPE(t_fftGrid) :: state, stateB, stateq, stateBq, StateDeriv, ekinGrid, chargeDen, rhomatGrid(4), statemq, stateBmq
      TYPE(t_fftGrid) :: stepFct
      INTEGER nu, iv, ir, istr, i, j
      INTEGER idens, ndens, ispin, iGp, iG, gVec(3), gInd
      REAL q0, q0_11, q0_22, norm, xk(3), fac
      REAL s, stateRadius, stateFFTRadius, stateFFTExtendedRadius
      REAL, PARAMETER:: tol_3 = 1.0e-3
      LOGICAL forw, l_dfpt, l_minusq

      ! local arrays
      INTEGER, ALLOCATABLE :: stateIndices(:)
      INTEGER, ALLOCATABLE :: stateBIndices(:)
      INTEGER, ALLOCATABLE :: stateqIndices(:)
      INTEGER, ALLOCATABLE :: stateBqIndices(:)
      INTEGER, ALLOCATABLE :: statemqIndices(:)
      INTEGER, ALLOCATABLE :: stateBmqIndices(:)
      INTEGER, ALLOCATABLE :: fieldSphereIndices(:)
      REAL wtf(ne), wtf1(ne)
      COMPLEX, ALLOCATABLE :: cwk(:), ecwk(:)

      ! per-thread objects of the loop over states
      TYPE(t_fft) :: fftBack, fftForw
      INTEGER :: nthr, fftBackend, glen
      LOGICAL :: l_obj, l_omp, l_gamma, l_stateCharge
      REAL, ALLOCATABLE    :: accR(:,:), ekR(:)
      COMPLEX, ALLOCATABLE :: accC(:,:), cwkS(:), starCharges(:), tempState(:)
      COMPLEX, ALLOCATABLE :: pwintStars(:), pwintAllStars(:)

      ! subroutines
      REAL dznrm2, dnrm2

      !------->          ABBREVIATIONS
      !
      !     ne    : number of occupied states
      !     nv    : number of g-components in eigenstate
      !     cwk   : complex work array: charge density in g-space (as stars)
      !     den%pw : charge density stored as stars
      !     we    : weights for the BZ-integration for a particular k-point
      !     omtil : volume (slab) unit cell, between -.5*D_tilde and +.5*D_tilde
      !     nstr  : number of members (arms) of reciprocal lattice (g) vector
      !             of each star.
      !     ng3_fft: number of stars in the  charge density  FFT-box
      !     ng3   : number of 3 dim. stars in the charge density sphere defined
      !             by gmax

      CALL timestart("pwden")

      l_dfpt = PRESENT(q_dfpt)

      l_minusq = PRESENT(lapwmq)

      stateRadius = MAXVAL(ABS(lapw%rk(:,:)))
      stateRadius = stateRadius + SQRT(DOT_PRODUCT(kpts%bk(:,ikpt),kpts%bk(:,ikpt)))
      IF (noco%l_noco) stateRadius = stateRadius + SQRT(DOT_PRODUCT(nococonv%qss(:),nococonv%qss(:)))
      IF (l_dfpt) stateRadius = stateRadius + SQRT(DOT_PRODUCT(q_dfpt,q_dfpt))

      stateFFTRadius = 2.0*stateRadius
      stateFFTExtendedRadius = 2.0*stateRadius

      IF (noco%l_noco.OR.noco%l_soc) THEN
         IF (banddos%dos .OR. banddos%vacdos .OR. input%cdinf.OR.banddos%band) THEN
            stateFFTExtendedRadius = 3.0*stateRadius+0.1
            CALL stepFct%init(cell,sym,stateFFTExtendedRadius+0.001)
            CALL stepFct%putFieldOnGrid(stars, stars%ustep, cell,stateFFTRadius+0.0005)
            CALL stepFct%fillFieldSphereIndexArray(stars, stateFFTRadius+0.0008, fieldSphereIndices)
            CALL priv_fft_dev(stepFct%dimensions, stepFct%grid, .FALSE., fieldSphereIndices)
         END IF
      END IF

      ALLOCATE (cwk(stars%ng3), ecwk(stars%ng3))

      ! Grids for the summed densities. All grids have the same dimensions, so they also provide
      ! the index arrays; the grids for the single states are set up per thread further below.
      IF (noco%l_noco) THEN
         DO i = 1, 4
            CALL rhomatGrid(i)%init(cell,sym,stateFFTExtendedRadius+0.001)
            rhomatGrid(i)%grid(:) = CMPLX(0.0,0.0)
         END DO
         CALL rhomatGrid(1)%fillFieldSphereIndexArray(stars, stateFFTRadius+0.0008, fieldSphereIndices)
         IF (noco%l_ss) THEN
            ALLOCATE(stateIndices(lapw%nv(1)))
            ALLOCATE(stateBIndices(lapw%nv(2)))
            CALL rhomatGrid(1)%fillStateIndexArray(lapw,1,stateIndices)
            CALL rhomatGrid(1)%fillStateIndexArray(lapw,2,stateBIndices)
            IF (l_dfpt) THEN
               ALLOCATE(stateqIndices(lapwq%nv(1)))
               ALLOCATE(stateBqIndices(lapwq%nv(2)))
               CALL rhomatGrid(1)%fillStateIndexArray(lapwq,1,stateqIndices)
               CALL rhomatGrid(1)%fillStateIndexArray(lapwq,2,stateBqIndices)
               IF (l_minusq) THEN
                  ALLOCATE(statemqIndices(lapwmq%nv(1)))
                  ALLOCATE(stateBmqIndices(lapwmq%nv(2)))
                  CALL rhomatGrid(1)%fillStateIndexArray(lapwmq,1,statemqIndices)
                  CALL rhomatGrid(1)%fillStateIndexArray(lapwmq,2,stateBmqIndices)
               END IF
            END IF
         ELSE
            ALLOCATE(stateIndices(lapw%nv(jspin)))
            ALLOCATE(stateBIndices(lapw%nv(jspin)))
            CALL rhomatGrid(1)%fillStateIndexArray(lapw,jspin,stateIndices)
            CALL rhomatGrid(1)%fillStateIndexArray(lapw,jspin,stateBIndices)
            IF (l_dfpt) THEN
               ALLOCATE(stateqIndices(lapwq%nv(jspin)))
               ALLOCATE(stateBqIndices(lapwq%nv(jspin)))
               CALL rhomatGrid(1)%fillStateIndexArray(lapwq,jspin,stateqIndices)
               CALL rhomatGrid(1)%fillStateIndexArray(lapwq,jspin,stateBqIndices)
               IF (l_minusq) THEN
                  ALLOCATE(statemqIndices(lapwmq%nv(jspin)))
                  ALLOCATE(stateBmqIndices(lapwmq%nv(jspin)))
                  CALL rhomatGrid(1)%fillStateIndexArray(lapwmq,jspin,statemqIndices)
                  CALL rhomatGrid(1)%fillStateIndexArray(lapwmq,jspin,stateBmqIndices)
               END IF
            END IF
         ENDIF
         glen = rhomatGrid(1)%gridLength
      ELSE
         CALL chargeDen%init(cell,sym,stateFFTExtendedRadius+0.001)
         chargeDen%grid(:) = CMPLX(0.0,0.0)
         CALL chargeDen%fillFieldSphereIndexArray(stars, stateFFTRadius+0.0008, fieldSphereIndices)
         IF (input%l_f) THEN
            CALL ekinGrid%init(cell,sym,stateFFTExtendedRadius+0.001)
            ekinGrid%grid(:) = CMPLX(0.0,0.0)
         END IF
         ALLOCATE(stateIndices(lapw%nv(jspin)))
         CALL chargeDen%fillStateIndexArray(lapw,jspin,stateIndices)
         IF (l_dfpt) THEN
            ALLOCATE(stateqIndices(lapwq%nv(jspin)))
            CALL chargeDen%fillStateIndexArray(lapwq,jspin,stateqIndices)
            IF (l_minusq) THEN
               ALLOCATE(statemqIndices(lapwmq%nv(jspin)))
               CALL chargeDen%fillStateIndexArray(lapwmq,jspin,statemqIndices)
            END IF
         END IF
         glen = chargeDen%gridLength
      ENDIF

      IF (.NOT.l_dfpt) THEN
         ! g=0 star: calculate the charge for this k-point and spin
         !           analytically to test the quality of FFT
         q0 = 0.0
         q0_11 = 0.0
         q0_22 = 0.0
         IF (noco%l_noco) THEN
            IF (.NOT. zmat%l_real) THEN
               DO nu = 1, ne
                  norm = dznrm2(lapw%nv(1),zMat%data_c(1:, nu),1) ! dznrm2 returns the 2-norm of the vector.
                  q0_11 = q0_11 + we(nu)*norm*norm
                  norm = dznrm2(lapw%nv(2),zMat%data_c(lapw%nv(1) + atoms%nlotot + 1:, nu),1) ! dznrm2 returns the 2-norm of the vector.
                  q0_22 = q0_22 + we(nu)*norm*norm
               ENDDO
            ENDIF
            q0_11 = q0_11/cell%omtil
            q0_22 = q0_22/cell%omtil
         ELSE
            IF (zmat%l_real) THEN
               DO nu = 1, ne
                  norm = dnrm2(lapw%nv(jspin),zMat%data_r(:, nu),1) ! dnrm2 returns the 2-norm of the vector.
                  q0 = q0 + we(nu)*norm*norm
               ENDDO
            ELSE
               DO nu = 1, ne
                  norm = dznrm2(lapw%nv(jspin),zMat%data_c(:, nu),1) ! dznrm2 returns the 2-norm of the vector.
                  q0 = q0 + we(nu)*norm*norm
               ENDDO
            ENDIF
            q0 = q0/cell%omtil
         ENDIF
      END IF

      wtf(:ne) = we(:ne)/cell%omtil
      IF (l_dfpt) wtf1(:ne) = we1(:ne)/cell%omtil
      fac = MERGE(1.0, 2.0, l_minusq) ! K: -q term = image of the +q term
      l_gamma = .FALSE.
      IF (l_dfpt) l_gamma = norm2(q_dfpt)<1e-8
      ndens = MERGE(4, 1, noco%l_noco)

      ! The charge of each state in the interstitial needs the step function integrals of the stars,
      ! which do not depend on the state
      l_stateCharge = banddos%dos .OR. banddos%vacdos .OR. input%cdinf .OR. banddos%band
      IF (l_stateCharge .AND. noco%l_noco) THEN
         ALLOCATE(pwintStars(stars%ng3_fft))
         DO istr = 1, stars%ng3_fft
            CALL pwint(stars, atoms, sym, cell, istr, pwintStars(istr))
         END DO
      END IF
      IF (l_stateCharge .AND. .NOT.noco%l_noco .AND. noco%l_soc .AND. input%jspins.EQ.2) THEN
         ALLOCATE(pwintAllStars(stars%ng3))
         pwintAllStars = CMPLX(0.0,0.0)
         CALL pwint_all(stars,atoms,sym ,cell,1,stars%ng3,pwintAllStars)
      END IF

#ifdef _OPENACC
      ! On the GPU the states are treated in batches with all data kept on the device
      CALL priv_states_gpu(stars, cell, input, noco, lapw, jspin, atoms%nlotot, ne, ev_list, ikpt, eig, wtf, wtf1, &
                           fac, l_dfpt, l_minusq, l_gamma, l_stateCharge, stateFFTRadius, zMat, &
                           stateIndices, stateBIndices, stateqIndices, stateBqIndices, statemqIndices, stateBmqIndices, &
                           stepFct, pwintStars, pwintAllStars, chargeDen, rhomatGrid, ekinGrid, dos, &
                           lapwq, zMat1, lapwmq, zMat1m)
#else
      ! The states are distributed over the OpenMP threads. Each thread has its own state grids,
      ! FFT plans and density accumulators, which are summed up at the end. This needs a thread-safe
      ! FFT: with the built-in FFT the loop runs on one thread, with SpFFT (plans depend on the
      ! index sets) the FFTs are done by fft_interface.
      fftBackend = selecFFT(.TRUE.)
      l_obj = fftBackend /= spFFT_const
      l_omp = (fftBackend == fftw_const) .OR. (fftBackend == mklFFT_const)
      nthr = 1
      !$ IF (l_omp) nthr = MAX(1, MIN(omp_get_max_threads(), ne))

      !$OMP PARALLEL NUM_THREADS(nthr) DEFAULT(shared) &
      !$OMP PRIVATE(nu,ir,istr,iv,i,j,s,xk,forw,state,stateB,stateq,stateBq,statemq,stateBmq,stateDeriv) &
      !$OMP PRIVATE(fftBack,fftForw,accR,accC,ekR,cwkS,starCharges,tempState)
      IF (noco%l_noco) THEN
         CALL state%init(cell,sym,stateFFTExtendedRadius+0.001)
         CALL stateB%init(cell,sym,stateFFTExtendedRadius+0.001)
         IF (l_dfpt) THEN
            CALL stateq%init(cell,sym,stateFFTExtendedRadius+0.001)
            CALL stateBq%init(cell,sym,stateFFTExtendedRadius+0.001)
            IF (l_minusq) THEN
               CALL statemq%init(cell,sym,stateFFTExtendedRadius+0.001)
               CALL stateBmq%init(cell,sym,stateFFTExtendedRadius+0.001)
            END IF
         END IF
      ELSE
         CALL state%init(cell,sym,stateFFTExtendedRadius+0.001)
         IF (l_dfpt) THEN
            CALL stateq%init(cell,sym,stateFFTExtendedRadius+0.001)
            IF (l_minusq) CALL statemq%init(cell,sym,stateFFTExtendedRadius+0.001)
         END IF
         IF (input%l_f) THEN
            CALL stateDeriv%init(cell,sym,stateFFTExtendedRadius+0.001)
            ALLOCATE(ekR(0:glen-1), tempState(lapw%nv(jspin)))
            ekR = 0.0
         END IF
      END IF
      ! Densities of the SCF case are real, the DFPT ones complex
      IF (l_dfpt) THEN
         ALLOCATE(accC(0:glen-1,ndens))
         accC = CMPLX(0.0,0.0)
      ELSE
         ALLOCATE(accR(0:glen-1,ndens))
         accR = 0.0
      END IF
      IF (l_stateCharge) ALLOCATE(cwkS(stars%ng3), starCharges(stars%ng3))
      IF (l_obj) THEN
         CALL fftBack%init(state%dimensions, .FALSE., stateIndices)
         CALL fftForw%init(state%dimensions, .TRUE., fieldSphereIndices)
      END IF

      ! LOOP OVER OCCUPIED STATES
      !$OMP DO SCHEDULE(dynamic)
      DO nu = 1, ne

         ! Do inverse FFT of state and add related charge density to overall charge density on real-space mesh.
         IF (noco%l_noco) THEN
            forw = .FALSE.
            IF (noco%l_ss) THEN
               CALL state%putComplexStateOnGrid(lapw, 1, zMat%data_c(1:lapw%nv(1),nu))
               CALL stateB%putComplexStateOnGrid(lapw, 2, zMat%data_c(lapw%nv(1) + atoms%nlotot+1:lapw%nv(1) + atoms%nlotot+lapw%nv(2),nu))
               IF (l_dfpt) THEN
                  CALL stateq%putComplexStateOnGrid(lapwq, 1, zMat1%data_c(1:lapwq%nv(1),nu))
                  CALL stateBq%putComplexStateOnGrid(lapwq, 2, zMat1%data_c(lapwq%nv(1) + atoms%nlotot+1:lapwq%nv(1) + atoms%nlotot+lapwq%nv(2),nu))
                  IF (l_minusq) THEN
                     CALL statemq%putComplexStateOnGrid(lapwmq, 1, zMat1m%data_c(1:lapwmq%nv(1),nu))
                     CALL stateBmq%putComplexStateOnGrid(lapwmq, 2, zMat1m%data_c(lapwmq%nv(1) + atoms%nlotot+1:lapwmq%nv(1) + atoms%nlotot+lapwmq%nv(2),nu))
                  END IF
               END IF
            ELSE
               CALL state%putComplexStateOnGrid(lapw, jspin, zMat%data_c(1:lapw%nv(jspin),nu))
               CALL stateB%putComplexStateOnGrid(lapw, jspin, zMat%data_c(lapw%nv(1) + atoms%nlotot+1:lapw%nv(1) + atoms%nlotot+lapw%nv(jspin),nu))
               IF (l_dfpt) THEN
                  CALL stateq%putComplexStateOnGrid(lapwq, jspin, zMat1%data_c(1:lapwq%nv(1),nu))
                  CALL stateBq%putComplexStateOnGrid(lapwq, jspin, zMat1%data_c(lapwq%nv(1) + atoms%nlotot+1:lapwq%nv(1) + atoms%nlotot+lapwq%nv(jspin),nu))
                  IF (l_minusq) THEN
                     CALL statemq%putComplexStateOnGrid(lapwmq, jspin, zMat1m%data_c(1:lapwmq%nv(1),nu))
                     CALL stateBmq%putComplexStateOnGrid(lapwmq, jspin, zMat1m%data_c(lapwmq%nv(1) + atoms%nlotot+1:lapwmq%nv(1) + atoms%nlotot+lapwmq%nv(jspin),nu))
                  END IF
               END IF
            END IF

            CALL priv_fft(fftBack, l_obj, state, forw, stateIndices)
            CALL priv_fft(fftBack, l_obj, stateB, forw, stateBIndices)
            IF (l_dfpt) THEN
               CALL priv_fft(fftBack, l_obj, stateq, forw, stateqIndices)
               CALL priv_fft(fftBack, l_obj, stateBq, forw, stateBqIndices)
               IF (l_minusq) THEN
                  CALL priv_fft(fftBack, l_obj, statemq, forw, statemqIndices)
                  CALL priv_fft(fftBack, l_obj, stateBmq, forw, stateBmqIndices)
               END IF
            END IF

            ! In the non-collinear case rho becomes a hermitian 2x2
            ! matrix (rhomatGrid).
            IF (.NOT.l_dfpt) THEN
            DO ir = 0, glen - 1
               !In this order: rho_11, rho_22, m_x/2, m_y/2
               accR(ir,1) = accR(ir,1) + wtf(nu) * (REAL(state%grid(ir))**2 + AIMAG(state%grid(ir))**2)
               accR(ir,2) = accR(ir,2) + wtf(nu) * (REAL(stateB%grid(ir))**2 + AIMAG(stateB%grid(ir))**2)
               accR(ir,3) = accR(ir,3) + wtf(nu) * (REAL(state%grid(ir))*REAL(stateB%grid(ir)) + AIMAG(state%grid(ir))*AIMAG(stateB%grid(ir)))
               accR(ir,4) = accR(ir,4) + wtf(nu) * (REAL(state%grid(ir))*AIMAG(stateB%grid(ir)) - AIMAG(state%grid(ir))*REAL(stateB%grid(ir)))
            END DO
            ELSE
               !DFPT NOCO
               DO ir = 0, glen - 1
                  !In this order: rho^1_11, rho^1_22, rho^1_12 (block (2,1)), rho^1_21 (block (1,2))
                  accC(ir,1) = accC(ir,1) + wtf(nu) * fac * CONJG(state%grid(ir))*stateq%grid(ir)
                  accC(ir,2) = accC(ir,2) + wtf(nu) * fac * CONJG(stateB%grid(ir))*stateBq%grid(ir)
                  accC(ir,3) = accC(ir,3) + wtf(nu) * fac * CONJG(stateB%grid(ir))*stateq%grid(ir)
                  accC(ir,4) = accC(ir,4) + wtf(nu) * fac * CONJG(state%grid(ir))*stateBq%grid(ir)
                  IF (l_minusq) THEN
                     accC(ir,1) = accC(ir,1) + wtf(nu) * CONJG(statemq%grid(ir))*state%grid(ir)
                     accC(ir,2) = accC(ir,2) + wtf(nu) * CONJG(stateBmq%grid(ir))*stateB%grid(ir)
                     accC(ir,3) = accC(ir,3) + wtf(nu) * CONJG(stateBmq%grid(ir))*state%grid(ir)
                     accC(ir,4) = accC(ir,4) + wtf(nu) * CONJG(statemq%grid(ir))*stateB%grid(ir)
                  END IF
                  IF (l_gamma) THEN
                     accC(ir,1) = accC(ir,1) + wtf1(nu) * ABS(state%grid(ir))**2
                     accC(ir,2) = accC(ir,2) + wtf1(nu) * ABS(stateB%grid(ir))**2
                     accC(ir,3) = accC(ir,3) + wtf1(nu) * CONJG(stateB%grid(ir))*state%grid(ir)
                     accC(ir,4) = accC(ir,4) + wtf1(nu) * CONJG(state%grid(ir))*stateB%grid(ir)
                  END IF
               END DO
            END IF

            ! In a non-collinear calculation the interstitial charge
            ! cannot be calculated by a simple substraction if the
            ! muffin-tin (and vacuum) charge is know, because the
            ! total charge does not need to be one in each spin-
            ! channel. Thus it has to be calculated explicitly, if
            ! it is needed.
            IF (l_stateCharge) THEN
               DO ir = 0, state%gridLength - 1
                  state%grid(ir) = ABS(state%grid(ir))**2
                  stateB%grid(ir) = ABS(stateB%grid(ir))**2
               END DO

               forw = .TRUE.
               CALL priv_fft(fftForw, l_obj, state, forw, fieldSphereIndices)
               CALL priv_fft(fftForw, l_obj, stateB, forw, fieldSphereIndices)

               CALL state%takeFieldFromGrid(stars, cwkS, stateFFTRadius+0.0005)
               DO istr = 1, stars%ng3
                  cwkS(istr) = REAL(stars%nstr(istr))*cwkS(istr)
               END DO
               DO istr = 1, stars%ng3_fft
                  dos%qis(ev_list(nu), ikpt, 1) = dos%qis(ev_list(nu), ikpt, 1) + REAL(cwkS(istr)*pwintStars(istr))/cell%omtil
                  dos%qTot(ev_list(nu), ikpt, 1) = dos%qTot(ev_list(nu), ikpt, 1) + REAL(cwkS(istr)*pwintStars(istr))/cell%omtil
               ENDDO

               CALL stateB%takeFieldFromGrid(stars, cwkS, stateFFTRadius+0.0005)
               DO istr = 1, stars%ng3
                  cwkS(istr) = REAL(stars%nstr(istr))*cwkS(istr)
               END DO
               DO istr = 1, stars%ng3_fft
                  dos%qis(ev_list(nu), ikpt, input%jspins) = dos%qis(ev_list(nu), ikpt, input%jspins) + REAL(cwkS(istr)*pwintStars(istr))/cell%omtil
                  dos%qTot(ev_list(nu), ikpt, input%jspins) = dos%qTot(ev_list(nu), ikpt, input%jspins) + REAL(cwkS(istr)*pwintStars(istr))/cell%omtil
               ENDDO
            ENDIF

         ELSE
            CALL state%putStateOnGrid(lapw, jSpin, zMat, nu)
            IF (l_dfpt) THEN
               CALL stateq%putStateOnGrid(lapwq, jspin, zMat1, nu)
               IF (l_minusq) THEN
                  CALL statemq%putStateOnGrid(lapwmq, jspin, zMat1m, nu)
               END IF
            END IF

            forw = .FALSE.
            ! The following FFT is a general complex to complex FFT
            ! For zmat%l_real this should be turned into a real to real FFT at some point
            ! Note: For the moment also no zero-indices for SpFFT provided
            CALL priv_fft(fftBack, l_obj, state, forw, stateIndices)
            IF (l_dfpt) THEN
               CALL priv_fft(fftBack, l_obj, stateq, forw, stateqIndices)
               IF (l_minusq) THEN
                  CALL priv_fft(fftBack, l_obj, statemq, forw, statemqIndices)
               END IF
            END IF
            IF (.NOT.l_dfpt) THEN
               DO ir = 0, glen - 1
                  accR(ir,1) = accR(ir,1) + wtf(nu) * (REAL(state%grid(ir))**2 + AIMAG(state%grid(ir))**2)
               END DO
            ELSE
               IF (.NOT.l_minusq) THEN
                  DO ir = 0, glen - 1
                     accC(ir,1) = accC(ir,1) + wtf(nu) * 2 * CONJG(state%grid(ir)) * stateq%grid(ir)
                     IF (l_gamma) accC(ir,1) = accC(ir,1) + wtf1(nu) * ABS(state%grid(ir))**2
                  END DO
               ELSE
                  DO ir = 0, glen - 1
                     accC(ir,1) = accC(ir,1) + wtf(nu) * (CONJG(statemq%grid(ir)) * state%grid(ir) &
                                                        + CONJG(state%grid(ir)) * stateq%grid(ir))
                     IF (l_gamma) accC(ir,1) = accC(ir,1) + wtf1(nu) * ABS(state%grid(ir))**2
                  END DO
               END IF
            END IF

            IF (input%l_f) THEN
               DO ir = 0, glen - 1
                  ekR(ir) = ekR(ir) - wtf(nu) * eig(nu) * (REAL(state%grid(ir))**2 + AIMAG(state%grid(ir))**2)
               END DO

               DO j = 1, 3
                  DO iv = 1, lapw%nv(jspin)
                     xk = lapw%gvec(:, iv, jspin) + lapw%bkpt
                     s = 0.0
                     DO i = 1, 3
                        s = s + xk(i)*cell%bmat(i, j)
                     ENDDO
                     IF (zmat%l_real) THEN
                        tempState(iv) = s*zMat%data_r(iv, nu)
                     ELSE
                        tempState(iv) = s*zMat%data_c(iv, nu)
                     END IF
                  END DO
                  CALL stateDeriv%putComplexStateOnGrid(lapw, jSpin, tempState)
                  CALL priv_fft(fftBack, l_obj, stateDeriv, forw, stateIndices)
                  DO ir = 0, glen - 1
                     ekR(ir) = ekR(ir) + wtf(nu) * 0.5 * (REAL(stateDeriv%grid(ir))**2 + AIMAG(stateDeriv%grid(ir))**2)
                  END DO
               END DO
            END IF

            IF (noco%l_soc.AND.input%jspins.EQ.2) THEN
               IF (l_stateCharge) THEN
                  DO ir = 0, state%gridLength - 1
                     state%grid(ir) = conjg(state%grid(ir)) * stepFct%grid(ir) * state%grid(ir)
                  END DO

                  forw = .TRUE.
                  CALL priv_fft(fftForw, l_obj, state, forw, fieldSphereIndices)

                  cwkS = CMPLX(0.0,0.0)
                  CALL state%takeFieldFromGrid(stars, cwkS, stateFFTRadius+0.0005)
                  starCharges(:) = pwintAllStars(:) * cwkS(:) * stars%nstr(:) / cell%omtil
                  dos%qis(ev_list(nu), ikpt, jSpin) = dos%qis(ev_list(nu), ikpt, jSpin) + REAL(SUM(starCharges(:)))
                  dos%qTot(ev_list(nu), ikpt, jSpin) = dos%qTot(ev_list(nu), ikpt, jSpin) + REAL(SUM(starCharges(:)))
               END IF
            END IF
         END IF
      END DO
      !$OMP END DO

      ! Sum up the densities of all threads
      !$OMP CRITICAL (pwden_sum)
      IF (noco%l_noco) THEN
         DO i = 1, 4
            IF (l_dfpt) THEN
               rhomatGrid(i)%grid(:) = rhomatGrid(i)%grid(:) + accC(:,i)
            ELSE
               rhomatGrid(i)%grid(:) = rhomatGrid(i)%grid(:) + accR(:,i)
            END IF
         END DO
      ELSE
         IF (l_dfpt) THEN
            chargeDen%grid(:) = chargeDen%grid(:) + accC(:,1)
         ELSE
            chargeDen%grid(:) = chargeDen%grid(:) + accR(:,1)
         END IF
         IF (input%l_f) ekinGrid%grid(:) = ekinGrid%grid(:) + ekR(:)
      END IF
      !$OMP END CRITICAL (pwden_sum)
      IF (l_obj) THEN
         CALL fftBack%free()
         CALL fftForw%free()
      END IF
      if (allocated(state%grid)) deallocate(state%grid)
      if (allocated(stateB%grid)) deallocate(stateB%grid)
      if (allocated(stateq%grid)) deallocate(stateq%grid)
      if (allocated(stateBq%grid)) deallocate(stateBq%grid)
      if (allocated(statemq%grid)) deallocate(statemq%grid)
      if (allocated(stateBmq%grid)) deallocate(stateBmq%grid)
      if (allocated(stateDeriv%grid)) deallocate(stateDeriv%grid)
      if (allocated(accR)) deallocate(accR)
      if (allocated(accC)) deallocate(accC)
      if (allocated(ekR)) deallocate(ekR)
      if (allocated(cwkS)) deallocate(cwkS)
      if (allocated(starCharges)) deallocate(starCharges)
      if (allocated(tempState)) deallocate(tempState)
      !$OMP END PARALLEL
#endif

      ! END OUTER LOOP OVER STATES NU

      ! Perform FFT transform of charge density to reciprocal space.

      ! In a collinear calculation pwden is called once per spin.
      ! However in a non-collinear calculation pwden is only called once
      ! and all four components of the density matrix (n_11 n_22 n_12
      ! n_21) have to be calculated at once.
      ndens = 1
      IF (noco%l_noco) ndens = 4
      DO idens = 1, ndens
         forw = .TRUE.
         IF (noco%l_noco) THEN
            CALL priv_fft_dev(rhomatGrid(idens)%dimensions, rhomatGrid(idens)%grid, forw, fieldSphereIndices)
         ELSE
            ! The following FFT is a general complex to complex FFT
            ! For zmat%l_real this should be turned into a real to real FFT at some point
            ! Note: For the moment also no zero-indices for SpFFT provided
            CALL priv_fft_dev(chargeDen%dimensions, chargeDen%grid, forw, fieldSphereIndices)
            IF (input%l_f) THEN
               CALL priv_fft_dev(ekinGrid%dimensions, ekinGrid%grid, forw, fieldSphereIndices)
            END IF
         ENDIF

         ! collect stars from the fft-grid
         cwk = 0.0
         ecwk = 0.0
         IF (noco%l_noco) THEN
            CALL rhomatGrid(idens)%takeFieldFromGrid(stars, cwk, stateFFTRadius+0.0005)
         ELSE
            CALL chargeDen%takeFieldFromGrid(stars, cwk, stateFFTRadius+0.0005)
            IF (input%l_f) THEN
               CALL ekinGrid%takeFieldFromGrid(stars, ecwk, stateFFTRadius+0.0005)
            END IF
         ENDIF

         IF (input%l_useapw) THEN
            IF (input%l_f) THEN
               CALL force_b8(atoms, ecwk, stars, sym, cell, jspin, results%force, f_b8)
            ENDIF
         ENDIF

         IF (.NOT.l_dfpt) THEN
         ! check charge neutralilty
         IF ((idens .EQ. 1) .OR. (idens .EQ. 2)) THEN
            IF (noco%l_noco) THEN
               IF (idens .EQ. 1) THEN
                  q0 = q0_11
               ELSE
                  q0 = q0_22
               ENDIF
            ENDIF
            IF (ABS(q0) .GT. 1.0e-9) THEN
               IF (ABS(q0 - REAL(cwk(1)))/q0 .GT. tol_3) THEN
                  WRITE (99, *) "XX:", ne, lapw%nv
                  IF (zmat%l_real) THEN
                     DO istr = 1, SIZE(zMat%data_r, 2)
                        WRITE (99, *) "X:", istr, zMat%data_r(:, istr)
                     ENDDO
                  ELSE
                     DO istr = 1, SIZE(zMat%data_c, 2)
                        WRITE (99, *) "X:", istr, zMat%data_c(:, istr)
                     ENDDO
                  ENDIF
                  WRITE (oUnit, '(''bad quality of charge density'',2f13.8)') q0, REAL(cwk(1))
                  CALL juDFT_warn('pwden: bad quality of charge')
               ENDIF
            ENDIF
         ENDIF
         END IF

         ! add charge density to existing one
         IF (idens .LE. 2) THEN
            ! add to spin-up or -down density (collinear & non-collinear)
            ispin = jspin
            IF (noco%l_noco) ispin = idens
            DO istr = 1, stars%ng3_fft
               den%pw(istr, ispin) = den%pw(istr, ispin) + cwk(istr)
            ENDDO
         ELSE IF (l_dfpt) THEN
            ! idens 3: block (2,1) = rho_12, idens 4: block (1,2) = rho_21
            IF (l_minusq) THEN
               DO istr = 1, stars%ng3_fft
                  den%pw(istr, idens) = den%pw(istr, idens) + cwk(istr)
               ENDDO
            ELSE
               ! K: both off-diagonal blocks are averaged into rho_12 and rho_21
               DO istr = 1, stars%ng3_fft
                  den%pw(istr, 3:4) = den%pw(istr, 3:4) + 0.5*cwk(istr)
               ENDDO
            END IF
         ELSE IF (idens .EQ. 3) THEN
            ! add to off-diag. part of density matrix (only non-collinear)
            DO istr = 1, stars%ng3_fft
               den%pw(istr, 3) = den%pw(istr, 3) + cwk(istr)
            ENDDO
            !IF (l_dfpt) THEN
            !   ! add to other off-diag. part of density matrix (only non-collinear)
            !   DO istr = 1, stars%ng3_fft
            !      den%pw(istr, 4) = den%pw(istr, 4) + cwk(istr)
            !   ENDDO
            !END IF
         ELSE
            ! add to off-diag. part of density matrix (only non-collinear)
            DO istr = 1, stars%ng3_fft
               den%pw(istr, 3) = den%pw(istr, 3) - ImagUnit*cwk(istr)
            ENDDO
            !IF (l_dfpt) THEN
            !   DO istr = 1, stars%ng3_fft
            !      den%pw(istr, 4) = den%pw(istr, 4) + ImagUnit*cwk(istr)
            !   ENDDO
            !END IF
         ENDIF
      ENDDO

      DEALLOCATE (cwk, ecwk)

      CALL timestop("pwden")

   END SUBROUTINE pwden

   SUBROUTINE priv_fft(fft, l_obj, fftGrid, forw, indices)
      !! FFT of a grid with a prepared plan (l_obj) or via fft_interface
      TYPE(t_fft),     INTENT(INOUT) :: fft
      LOGICAL,         INTENT(IN)    :: l_obj, forw
      TYPE(t_fftGrid), INTENT(INOUT) :: fftGrid
      INTEGER,         INTENT(IN)    :: indices(:)

      IF (l_obj) THEN
         CALL timestart("FFT_interface")
         CALL fft%exec(fftGrid%grid)
         CALL timestop("FFT_interface")
      ELSE
         CALL fft_interface(3, fftGrid%dimensions(:), fftGrid%grid, forw, indices)
      END IF
   END SUBROUTINE priv_fft

   SUBROUTINE priv_fft_dev(dims, grid, forw, indices)
      !! FFT of a single grid; in OpenACC builds the grid is put on the device, so that cuFFT is used
      INTEGER, INTENT(IN)    :: dims(3)
      COMPLEX, INTENT(INOUT) :: grid(:)
      LOGICAL, INTENT(IN)    :: forw
      INTEGER, INTENT(IN)    :: indices(:)

      !$acc data copy(grid)
      CALL fft_interface(3, dims, grid, forw, indices)
      !$acc end data
   END SUBROUTINE priv_fft_dev

   SUBROUTINE priv_states_gpu(stars, cell, input, noco, lapw, jspin, nlotot, ne, ev_list, ikpt, eig, wtf, wtf1, &
                              fac, l_dfpt, l_minusq, l_gamma, l_stateCharge, stateFFTRadius, zMat, &
                              stateIndices, stateBIndices, stateqIndices, stateBqIndices, statemqIndices, stateBmqIndices, &
                              stepFct, pwintStars, pwintAllStars, chargeDen, rhomatGrid, ekinGrid, dos, &
                              lapwq, zMat1, lapwmq, zMat1m)
      !! Densities of the states in the interstitial with the states treated in batches on the device:
      !! the coefficients of a batch are put on the FFT grids, transformed to real space by a batched FFT
      !! and accumulated into the densities. Only the summed densities and the charges of the states
      !! are transferred back.
      TYPE(t_stars),   INTENT(IN)    :: stars
      TYPE(t_cell),    INTENT(IN)    :: cell
      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_noco),    INTENT(IN)    :: noco
      TYPE(t_lapw),    INTENT(IN)    :: lapw
      INTEGER,         INTENT(IN)    :: jspin, nlotot, ne, ev_list(:), ikpt
      REAL,            INTENT(IN)    :: eig(:), wtf(:), wtf1(:), fac, stateFFTRadius
      LOGICAL,         INTENT(IN)    :: l_dfpt, l_minusq, l_gamma, l_stateCharge
      TYPE(t_mat),     INTENT(IN)    :: zMat
      INTEGER, ALLOCATABLE, INTENT(IN) :: stateIndices(:), stateBIndices(:), stateqIndices(:), stateBqIndices(:)
      INTEGER, ALLOCATABLE, INTENT(IN) :: statemqIndices(:), stateBmqIndices(:)
      TYPE(t_fftGrid), INTENT(IN)    :: stepFct
      COMPLEX, ALLOCATABLE, INTENT(IN) :: pwintStars(:), pwintAllStars(:)
      TYPE(t_fftGrid), INTENT(INOUT) :: chargeDen, rhomatGrid(4), ekinGrid
      TYPE(t_dos),     INTENT(INOUT) :: dos
      TYPE(t_lapw), OPTIONAL, INTENT(IN) :: lapwq, lapwmq
      TYPE(t_mat),  OPTIONAL, INTENT(IN) :: zMat1, zMat1m

      COMPLEX, ALLOCATABLE :: psiA(:,:), psiB(:,:), psiAq(:,:), psiBq(:,:), psiAmq(:,:), psiBmq(:,:), psiD(:,:)
      COMPLEX, ALLOCATABLE :: accC(:,:), wq(:), step(:)
      REAL,    ALLOCATABLE :: accR(:,:), ekR(:), scal(:,:), qs(:,:)
      INTEGER, ALLOCATABLE :: gq(:)
      TYPE(t_fft) :: fftBack, fftForw
      INTEGER :: dims(3), glen, ndens, ngrids, nb, nbb, nu0, nu, ib, ir, iv, i, j, nvmax, offB, offBq, offBmq
      REAL    :: w, w1, s1, s2, s3, s4
      COMPLEX :: a, b, aq, bq, amq, bmq, c1, c2, c3, c4
      LOGICAL :: l_noco, l_f, l_qNoco, l_qSoc

      CALL timestart("pwden: states on device")
      l_noco  = noco%l_noco
      l_f     = input%l_f .AND. .NOT.l_noco
      l_qNoco = l_stateCharge .AND. l_noco
      l_qSoc  = l_stateCharge .AND. .NOT.l_noco .AND. noco%l_soc .AND. input%jspins == 2
      ndens   = MERGE(4, 1, l_noco)
      IF (l_noco) THEN
         dims = rhomatGrid(1)%dimensions
         glen = rhomatGrid(1)%gridLength
      ELSE
         dims = chargeDen%dimensions
         glen = chargeDen%gridLength
      END IF

      ngrids = MERGE(2, 1, l_noco)
      IF (l_dfpt) ngrids = ngrids*MERGE(3, 2, l_minusq)
      IF (l_f) ngrids = ngrids + 1
      nb = priv_batch_size(ne, glen, ngrids)

      ! Row offsets of the second spin block in the noco eigenvectors
      offB = lapw%nv(1) + nlotot
      offBq = 0
      offBmq = 0
      IF (l_dfpt .AND. l_noco) offBq = lapwq%nv(1) + nlotot
      IF (l_minusq .AND. l_noco) offBmq = lapwmq%nv(1) + nlotot

      ! Scaling of the coefficients: column 0 for the states, 1-3 for the gradients (forces)
      nvmax = SIZE(stateIndices)
      IF (ALLOCATED(stateBIndices)) nvmax = MAX(nvmax, SIZE(stateBIndices))
      IF (ALLOCATED(stateqIndices)) nvmax = MAX(nvmax, SIZE(stateqIndices))
      IF (ALLOCATED(stateBqIndices)) nvmax = MAX(nvmax, SIZE(stateBqIndices))
      IF (ALLOCATED(statemqIndices)) nvmax = MAX(nvmax, SIZE(statemqIndices))
      IF (ALLOCATED(stateBmqIndices)) nvmax = MAX(nvmax, SIZE(stateBmqIndices))
      ALLOCATE(scal(nvmax, 0:3))
      scal = 1.0
      IF (l_f) THEN
         DO j = 1, 3
            DO iv = 1, lapw%nv(jspin)
               scal(iv, j) = DOT_PRODUCT(lapw%gvec(:, iv, jspin) + lapw%bkpt, cell%bmat(:, j))
            END DO
         END DO
      END IF

      ! Grid indices and weights for the charges of the single states
      IF (l_qNoco) CALL priv_star_weights(stars, cell, dims, stateFFTRadius+0.0005, pwintStars, stars%ng3_fft, gq, wq)
      IF (l_qSoc) THEN
         CALL priv_star_weights(stars, cell, dims, stateFFTRadius+0.0005, pwintAllStars, stars%ng3, gq, wq)
         ALLOCATE(step(0:glen-1))
         step(:) = stepFct%grid(:)
      END IF

      ! Device data
      CALL priv_z_device(zMat, .TRUE.)
      IF (l_dfpt) CALL priv_z_device(zMat1, .TRUE.)
      IF (l_minusq) CALL priv_z_device(zMat1m, .TRUE.)
      !$acc enter data copyin(stateIndices, eig, wtf, wtf1, scal)
      IF (ALLOCATED(stateBIndices)) THEN
         !$acc enter data copyin(stateBIndices)
      END IF
      IF (ALLOCATED(stateqIndices)) THEN
         !$acc enter data copyin(stateqIndices)
      END IF
      IF (ALLOCATED(stateBqIndices)) THEN
         !$acc enter data copyin(stateBqIndices)
      END IF
      IF (ALLOCATED(statemqIndices)) THEN
         !$acc enter data copyin(statemqIndices)
      END IF
      IF (ALLOCATED(stateBmqIndices)) THEN
         !$acc enter data copyin(stateBmqIndices)
      END IF

      ALLOCATE(psiA(0:glen-1, nb))
      !$acc enter data create(psiA)
      IF (l_noco) THEN
         ALLOCATE(psiB(0:glen-1, nb))
         !$acc enter data create(psiB)
      END IF
      IF (l_dfpt) THEN
         ALLOCATE(psiAq(0:glen-1, nb))
         !$acc enter data create(psiAq)
         IF (l_noco) THEN
            ALLOCATE(psiBq(0:glen-1, nb))
            !$acc enter data create(psiBq)
         END IF
         ! Without -q the accumulation kernels still name psiAmq/psiBmq: placeholders keep them present
         ALLOCATE(psiAmq(0:MERGE(glen-1, 0, l_minusq), MERGE(nb, 1, l_minusq)))
         !$acc enter data create(psiAmq)
         IF (l_noco) THEN
            ALLOCATE(psiBmq(0:MERGE(glen-1, 0, l_minusq), MERGE(nb, 1, l_minusq)))
            !$acc enter data create(psiBmq)
         END IF
         ALLOCATE(accC(0:glen-1, ndens))
         accC = CMPLX(0.0, 0.0)
         !$acc enter data copyin(accC)
      ELSE
         ALLOCATE(accR(0:glen-1, ndens))
         accR = 0.0
         !$acc enter data copyin(accR)
      END IF
      IF (l_f) THEN
         ALLOCATE(psiD(0:glen-1, nb), ekR(0:glen-1))
         ekR = 0.0
         !$acc enter data create(psiD) copyin(ekR)
      END IF
      IF (l_qNoco .OR. l_qSoc) THEN
         ALLOCATE(qs(nb, 2))
         !$acc enter data copyin(gq, wq) create(qs)
      END IF
      IF (l_qSoc) THEN
         !$acc enter data copyin(step)
      END IF

      CALL fftBack%init(dims, .FALSE., batch_size=nb, l_gpu=.TRUE.)
      IF (l_qNoco .OR. l_qSoc) CALL fftForw%init(dims, .TRUE., batch_size=nb, l_gpu=.TRUE.)

      DO nu0 = 1, ne, nb
         nbb = MIN(nb, ne - nu0 + 1)

         ! States of the batch in real space
         CALL priv_put_batch(psiA, stateIndices, zMat, 0, nu0, nbb, scal(:, 0))
         CALL fftBack%exec_batch(psiA)
         IF (l_noco) THEN
            CALL priv_put_batch(psiB, stateBIndices, zMat, offB, nu0, nbb, scal(:, 0))
            CALL fftBack%exec_batch(psiB)
         END IF
         IF (l_dfpt) THEN
            CALL priv_put_batch(psiAq, stateqIndices, zMat1, 0, nu0, nbb, scal(:, 0))
            CALL fftBack%exec_batch(psiAq)
            IF (l_noco) THEN
               CALL priv_put_batch(psiBq, stateBqIndices, zMat1, offBq, nu0, nbb, scal(:, 0))
               CALL fftBack%exec_batch(psiBq)
            END IF
            IF (l_minusq) THEN
               CALL priv_put_batch(psiAmq, statemqIndices, zMat1m, 0, nu0, nbb, scal(:, 0))
               CALL fftBack%exec_batch(psiAmq)
               IF (l_noco) THEN
                  CALL priv_put_batch(psiBmq, stateBmqIndices, zMat1m, offBmq, nu0, nbb, scal(:, 0))
                  CALL fftBack%exec_batch(psiBmq)
               END IF
            END IF
         END IF

         ! Accumulation of the densities, in the same form as in the loop over the states on the host
         IF (l_noco .AND. .NOT.l_dfpt) THEN
            ! rho_11, rho_22, m_x/2, m_y/2
            !$acc parallel loop default(present) private(s1, s2, s3, s4, w, a, b)
            DO ir = 0, glen - 1
               s1 = 0.0
               s2 = 0.0
               s3 = 0.0
               s4 = 0.0
               !$acc loop seq
               DO ib = 1, nbb
                  w = wtf(nu0+ib-1)
                  a = psiA(ir, ib)
                  b = psiB(ir, ib)
                  s1 = s1 + w*(REAL(a)**2 + AIMAG(a)**2)
                  s2 = s2 + w*(REAL(b)**2 + AIMAG(b)**2)
                  s3 = s3 + w*(REAL(a)*REAL(b) + AIMAG(a)*AIMAG(b))
                  s4 = s4 + w*(REAL(a)*AIMAG(b) - AIMAG(a)*REAL(b))
               END DO
               accR(ir, 1) = accR(ir, 1) + s1
               accR(ir, 2) = accR(ir, 2) + s2
               accR(ir, 3) = accR(ir, 3) + s3
               accR(ir, 4) = accR(ir, 4) + s4
            END DO
         ELSE IF (l_noco) THEN
            ! rho^1_11, rho^1_22, rho^1_12 (block (2,1)), rho^1_21 (block (1,2))
            !$acc parallel loop default(present) private(c1, c2, c3, c4, w, w1, a, b, aq, bq, amq, bmq)
            DO ir = 0, glen - 1
               c1 = CMPLX(0.0, 0.0)
               c2 = CMPLX(0.0, 0.0)
               c3 = CMPLX(0.0, 0.0)
               c4 = CMPLX(0.0, 0.0)
               !$acc loop seq
               DO ib = 1, nbb
                  w = wtf(nu0+ib-1)
                  w1 = wtf1(nu0+ib-1)
                  a = psiA(ir, ib)
                  b = psiB(ir, ib)
                  aq = psiAq(ir, ib)
                  bq = psiBq(ir, ib)
                  c1 = c1 + w*fac*CONJG(a)*aq
                  c2 = c2 + w*fac*CONJG(b)*bq
                  c3 = c3 + w*fac*CONJG(b)*aq
                  c4 = c4 + w*fac*CONJG(a)*bq
                  IF (l_minusq) THEN
                     amq = psiAmq(ir, ib)
                     bmq = psiBmq(ir, ib)
                     c1 = c1 + w*CONJG(amq)*a
                     c2 = c2 + w*CONJG(bmq)*b
                     c3 = c3 + w*CONJG(bmq)*a
                     c4 = c4 + w*CONJG(amq)*b
                  END IF
                  IF (l_gamma) THEN
                     c1 = c1 + w1*ABS(a)**2
                     c2 = c2 + w1*ABS(b)**2
                     c3 = c3 + w1*CONJG(b)*a
                     c4 = c4 + w1*CONJG(a)*b
                  END IF
               END DO
               accC(ir, 1) = accC(ir, 1) + c1
               accC(ir, 2) = accC(ir, 2) + c2
               accC(ir, 3) = accC(ir, 3) + c3
               accC(ir, 4) = accC(ir, 4) + c4
            END DO
         ELSE IF (.NOT.l_dfpt) THEN
            !$acc parallel loop default(present) private(s1, a)
            DO ir = 0, glen - 1
               s1 = 0.0
               !$acc loop seq
               DO ib = 1, nbb
                  a = psiA(ir, ib)
                  s1 = s1 + wtf(nu0+ib-1)*(REAL(a)**2 + AIMAG(a)**2)
               END DO
               accR(ir, 1) = accR(ir, 1) + s1
            END DO
         ELSE
            !$acc parallel loop default(present) private(c1, w, w1, a, aq, amq)
            DO ir = 0, glen - 1
               c1 = CMPLX(0.0, 0.0)
               !$acc loop seq
               DO ib = 1, nbb
                  w = wtf(nu0+ib-1)
                  w1 = wtf1(nu0+ib-1)
                  a = psiA(ir, ib)
                  aq = psiAq(ir, ib)
                  IF (l_minusq) THEN
                     amq = psiAmq(ir, ib)
                     c1 = c1 + w*(CONJG(amq)*a + CONJG(a)*aq)
                  ELSE
                     c1 = c1 + w*2*CONJG(a)*aq
                  END IF
                  IF (l_gamma) c1 = c1 + w1*ABS(a)**2
               END DO
               accC(ir, 1) = accC(ir, 1) + c1
            END DO
         END IF

         ! Kinetic energy density for the forces
         IF (l_f) THEN
            !$acc parallel loop default(present) private(s1, a)
            DO ir = 0, glen - 1
               s1 = 0.0
               !$acc loop seq
               DO ib = 1, nbb
                  a = psiA(ir, ib)
                  s1 = s1 - wtf(nu0+ib-1)*eig(nu0+ib-1)*(REAL(a)**2 + AIMAG(a)**2)
               END DO
               ekR(ir) = ekR(ir) + s1
            END DO
            DO j = 1, 3
               CALL priv_put_batch(psiD, stateIndices, zMat, 0, nu0, nbb, scal(:, j))
               CALL fftBack%exec_batch(psiD)
               !$acc parallel loop default(present) private(s1, a)
               DO ir = 0, glen - 1
                  s1 = 0.0
                  !$acc loop seq
                  DO ib = 1, nbb
                     a = psiD(ir, ib)
                     s1 = s1 + wtf(nu0+ib-1)*0.5*(REAL(a)**2 + AIMAG(a)**2)
                  END DO
                  ekR(ir) = ekR(ir) + s1
               END DO
            END DO
         END IF

         ! Interstitial charges of the single states
         IF (l_qNoco) THEN
            !$acc parallel loop collapse(2) default(present)
            DO ib = 1, nb
               DO ir = 0, glen - 1
                  psiA(ir, ib) = ABS(psiA(ir, ib))**2
                  psiB(ir, ib) = ABS(psiB(ir, ib))**2
               END DO
            END DO
            CALL fftForw%exec_batch(psiA)
            CALL fftForw%exec_batch(psiB)
            CALL priv_qcharge(psiA, gq, wq, nbb, qs(:, 1))
            CALL priv_qcharge(psiB, gq, wq, nbb, qs(:, 2))
            DO ib = 1, nbb
               nu = ev_list(nu0+ib-1)
               dos%qis(nu, ikpt, 1) = dos%qis(nu, ikpt, 1) + qs(ib, 1)
               dos%qTot(nu, ikpt, 1) = dos%qTot(nu, ikpt, 1) + qs(ib, 1)
               dos%qis(nu, ikpt, input%jspins) = dos%qis(nu, ikpt, input%jspins) + qs(ib, 2)
               dos%qTot(nu, ikpt, input%jspins) = dos%qTot(nu, ikpt, input%jspins) + qs(ib, 2)
            END DO
         ELSE IF (l_qSoc) THEN
            !$acc parallel loop collapse(2) default(present)
            DO ib = 1, nb
               DO ir = 0, glen - 1
                  psiA(ir, ib) = CONJG(psiA(ir, ib))*step(ir)*psiA(ir, ib)
               END DO
            END DO
            CALL fftForw%exec_batch(psiA)
            CALL priv_qcharge(psiA, gq, wq, nbb, qs(:, 1))
            DO ib = 1, nbb
               nu = ev_list(nu0+ib-1)
               dos%qis(nu, ikpt, jspin) = dos%qis(nu, ikpt, jspin) + qs(ib, 1)
               dos%qTot(nu, ikpt, jspin) = dos%qTot(nu, ikpt, jspin) + qs(ib, 1)
            END DO
         END IF
      END DO

      ! Summed densities back to the host
      IF (l_dfpt) THEN
         !$acc exit data copyout(accC)
         IF (l_noco) THEN
            DO i = 1, 4
               rhomatGrid(i)%grid(:) = rhomatGrid(i)%grid(:) + accC(:, i)
            END DO
         ELSE
            chargeDen%grid(:) = chargeDen%grid(:) + accC(:, 1)
         END IF
      ELSE
         !$acc exit data copyout(accR)
         IF (l_noco) THEN
            DO i = 1, 4
               rhomatGrid(i)%grid(:) = rhomatGrid(i)%grid(:) + accR(:, i)
            END DO
         ELSE
            chargeDen%grid(:) = chargeDen%grid(:) + accR(:, 1)
         END IF
      END IF
      IF (l_f) THEN
         !$acc exit data copyout(ekR) delete(psiD)
         ekinGrid%grid(:) = ekinGrid%grid(:) + ekR(:)
      END IF

      CALL fftBack%free()
      IF (l_qNoco .OR. l_qSoc) THEN
         CALL fftForw%free()
         !$acc exit data delete(gq, wq, qs)
      END IF
      IF (l_qSoc) THEN
         !$acc exit data delete(step)
      END IF
      IF (ALLOCATED(psiBmq)) THEN
         !$acc exit data delete(psiBmq)
      END IF
      IF (ALLOCATED(psiAmq)) THEN
         !$acc exit data delete(psiAmq)
      END IF
      IF (ALLOCATED(psiBq)) THEN
         !$acc exit data delete(psiBq)
      END IF
      IF (ALLOCATED(psiAq)) THEN
         !$acc exit data delete(psiAq)
      END IF
      IF (ALLOCATED(psiB)) THEN
         !$acc exit data delete(psiB)
      END IF
      !$acc exit data delete(psiA)
      IF (ALLOCATED(stateBmqIndices)) THEN
         !$acc exit data delete(stateBmqIndices)
      END IF
      IF (ALLOCATED(statemqIndices)) THEN
         !$acc exit data delete(statemqIndices)
      END IF
      IF (ALLOCATED(stateBqIndices)) THEN
         !$acc exit data delete(stateBqIndices)
      END IF
      IF (ALLOCATED(stateqIndices)) THEN
         !$acc exit data delete(stateqIndices)
      END IF
      IF (ALLOCATED(stateBIndices)) THEN
         !$acc exit data delete(stateBIndices)
      END IF
      !$acc exit data delete(stateIndices, eig, wtf, wtf1, scal)
      IF (l_minusq) CALL priv_z_device(zMat1m, .FALSE.)
      IF (l_dfpt) CALL priv_z_device(zMat1, .FALSE.)
      CALL priv_z_device(zMat, .FALSE.)
      CALL timestop("pwden: states on device")
   END SUBROUTINE priv_states_gpu

   INTEGER FUNCTION priv_batch_size(ne, glen, ngrids)
      !! Number of states per batch, limited by the free device memory
      INTEGER, INTENT(IN) :: ne, glen, ngrids
#ifdef _OPENACC
      INTEGER(8) :: freeMem

      freeMem = acc_get_property(0, acc_device_current, acc_property_free_memory)
      ! half of the free memory; one more grid per state for the cuFFT work space
      priv_batch_size = INT(0.5*REAL(freeMem)/(16.0*REAL(glen)*REAL(ngrids+1)))
#else
      priv_batch_size = 16
#endif
      priv_batch_size = MAX(1, MIN(priv_batch_size, 64, ne))
   END FUNCTION priv_batch_size

   SUBROUTINE priv_star_weights(stars, cell, dims, gCutoff, pw, nmax, gq, wq)
      !! Grid points and weights with REAL(SUM(grid(gq)*wq)) = interstitial charge of a density given on the
      !! grid in reciprocal space (unnormalized forward FFT), from the stars up to nmax and their integrals pw
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_cell),  INTENT(IN) :: cell
      INTEGER,       INTENT(IN) :: dims(3), nmax
      REAL,          INTENT(IN) :: gCutoff
      COMPLEX,       INTENT(IN) :: pw(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: gq(:)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: wq(:)

      INTEGER, ALLOCATABLE :: gTmp(:)
      COMPLEX, ALLOCATABLE :: wTmp(:)
      INTEGER :: x, y, z, iStar, n
      REAL    :: fct

      ALLOCATE(gTmp((2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)))
      ALLOCATE(wTmp(SIZE(gTmp)))
      fct = 1.0/(REAL(PRODUCT(dims))*cell%omtil)
      n = 0
      DO z = -stars%mx3, stars%mx3
         DO y = -stars%mx2, stars%mx2
            DO x = -stars%mx1, stars%mx1
               iStar = stars%ig(x, y, z)
               IF (iStar == 0 .OR. iStar > nmax) CYCLE
               IF (stars%sk3(iStar) > gCutoff) CYCLE
               n = n + 1
               gTmp(n) = MODULO(x, dims(1)) + dims(1)*MODULO(y, dims(2)) + dims(1)*dims(2)*MODULO(z, dims(3))
               wTmp(n) = CONJG(stars%rgphs(x, y, z))*pw(iStar)*fct
            END DO
         END DO
      END DO
      gq = gTmp(1:n)
      wq = wTmp(1:n)
   END SUBROUTINE priv_star_weights

   SUBROUTINE priv_qcharge(psi, gq, wq, nbb, q)
      !! Charges REAL(SUM(psi(gq,ib)*wq)) of the first nbb grids of psi (device data), returned on the host
      COMPLEX, INTENT(IN)    :: psi(0:,:), wq(:)
      INTEGER, INTENT(IN)    :: gq(:), nbb
      REAL,    INTENT(INOUT) :: q(:)

      INTEGER :: ib, ip, np
      REAL    :: s

      np = SIZE(gq)
      !$acc parallel loop gang default(present) private(s)
      DO ib = 1, nbb
         s = 0.0
         !$acc loop vector reduction(+:s)
         DO ip = 1, np
            s = s + REAL(psi(gq(ip), ib)*wq(ip))
         END DO
         q(ib) = s
      END DO
      !$acc update self(q(1:nbb))
   END SUBROUTINE priv_qcharge

   SUBROUTINE priv_z_device(zmat, l_enter)
      !! Eigenvector coefficients to the device (l_enter) or removed from it
      TYPE(t_mat), INTENT(IN) :: zmat
      LOGICAL,     INTENT(IN) :: l_enter

      IF (zmat%l_real) THEN
         CALL priv_r_device(zmat%data_r, l_enter)
      ELSE
         CALL priv_c_device(zmat%data_c, l_enter)
      END IF
   END SUBROUTINE priv_z_device

   SUBROUTINE priv_r_device(z, l_enter)
      REAL,    INTENT(IN) :: z(:,:)
      LOGICAL, INTENT(IN) :: l_enter

      IF (l_enter) THEN
         !$acc enter data copyin(z)
      ELSE
         !$acc exit data delete(z)
      END IF
   END SUBROUTINE priv_r_device

   SUBROUTINE priv_c_device(z, l_enter)
      COMPLEX, INTENT(IN) :: z(:,:)
      LOGICAL, INTENT(IN) :: l_enter

      IF (l_enter) THEN
         !$acc enter data copyin(z)
      ELSE
         !$acc exit data delete(z)
      END IF
   END SUBROUTINE priv_c_device

   SUBROUTINE priv_put_batch(psi, idx, zmat, off, nu0, nbb, scal)
      !! Coefficients off+1..off+SIZE(idx) of the states nu0..nu0+nbb-1, scaled by scal, on the FFT grids
      !! psi(:,1:nbb) (device data); all other grid points are zero
      COMPLEX,     INTENT(INOUT) :: psi(0:,:)
      INTEGER,     INTENT(IN)    :: idx(:), off, nu0, nbb
      TYPE(t_mat), INTENT(IN)    :: zmat
      REAL,        INTENT(IN)    :: scal(:)

      IF (zmat%l_real) THEN
         CALL priv_put_r(psi, idx, zmat%data_r, off, nu0, nbb, scal)
      ELSE
         CALL priv_put_c(psi, idx, zmat%data_c, off, nu0, nbb, scal)
      END IF
   END SUBROUTINE priv_put_batch

   SUBROUTINE priv_put_r(psi, idx, z, off, nu0, nbb, scal)
      COMPLEX, INTENT(INOUT) :: psi(0:,:)
      INTEGER, INTENT(IN)    :: idx(:), off, nu0, nbb
      REAL,    INTENT(IN)    :: z(:,:), scal(:)

      INTEGER :: ib, ir, iv, glen, nb, nv

      glen = SIZE(psi, 1)
      nb = SIZE(psi, 2)
      nv = SIZE(idx)
      !$acc parallel loop collapse(2) default(present)
      DO ib = 1, nb
         DO ir = 0, glen - 1
            psi(ir, ib) = CMPLX(0.0, 0.0)
         END DO
      END DO
      !$acc parallel loop collapse(2) default(present)
      DO ib = 1, nbb
         DO iv = 1, nv
            psi(idx(iv), ib) = CMPLX(z(off+iv, nu0+ib-1)*scal(iv), 0.0)
         END DO
      END DO
   END SUBROUTINE priv_put_r

   SUBROUTINE priv_put_c(psi, idx, z, off, nu0, nbb, scal)
      COMPLEX, INTENT(INOUT) :: psi(0:,:)
      INTEGER, INTENT(IN)    :: idx(:), off, nu0, nbb
      COMPLEX, INTENT(IN)    :: z(:,:)
      REAL,    INTENT(IN)    :: scal(:)

      INTEGER :: ib, ir, iv, glen, nb, nv

      glen = SIZE(psi, 1)
      nb = SIZE(psi, 2)
      nv = SIZE(idx)
      !$acc parallel loop collapse(2) default(present)
      DO ib = 1, nb
         DO ir = 0, glen - 1
            psi(ir, ib) = CMPLX(0.0, 0.0)
         END DO
      END DO
      !$acc parallel loop collapse(2) default(present)
      DO ib = 1, nbb
         DO iv = 1, nv
            psi(idx(iv), ib) = z(off+iv, nu0+ib-1)*scal(iv)
         END DO
      END DO
   END SUBROUTINE priv_put_c
END MODULE m_pwden
