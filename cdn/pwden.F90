!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_pwden
CONTAINS
   SUBROUTINE pwden(stars, kpts, banddos, oneD, input, fmpi, noco, nococonv, cell, atoms, sym, &
                    ikpt, jspin, lapw, ne, ev_list, we, eig, den, results, f_b8, zMat, dos)
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

!DEC$ NOOPTIMIZE
#include"cpp_double.h"
      USE m_types
      USE m_types_dos
      USE m_constants
      USE m_forceb8
      USE m_pwint
      USE m_juDFT
      USE m_types_fftGrid
      USE m_fft_interface

      IMPLICIT NONE

      TYPE(t_lapw), INTENT(IN)       :: lapw
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      TYPE(t_oneD), INTENT(IN)       :: oneD
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
      INTEGER, INTENT(IN)    :: ev_list(ne)
      INTEGER, INTENT(IN)    :: ne
      INTEGER, INTENT(IN)    :: ikpt
      INTEGER, INTENT(IN)    :: jspin
      COMPLEX, INTENT(INOUT) ::  f_b8(3, atoms%ntype)

      ! local variables
      TYPE(t_fftGrid) :: state, stateB, StateDeriv, ekinGrid, chargeDen, rhomatGrid(4)
      TYPE(t_fftGrid) :: stepFct
      INTEGER nu, iv, ir, istr, i, j
      INTEGER idens, ndens, ispin
      REAL q0, q0_11, q0_22, norm, xk(3)
      REAL s, stateRadius, stateFFTRadius, stateFFTExtendedRadius
      COMPLEX x
      REAL, PARAMETER:: tol_3 = 1.0e-3
      LOGICAL forw

      ! local arrays
      INTEGER, ALLOCATABLE :: stateIndices(:)
      INTEGER, ALLOCATABLE :: stateBIndices(:)
      INTEGER, ALLOCATABLE :: fieldSphereIndices(:)
      REAL wtf(ne)
      COMPLEX tempState(lapw%nv(jspin)), starCharges(stars%ng3)
      COMPLEX, ALLOCATABLE :: cwk(:), ecwk(:)

      ! subroutines
      REAL dznrm2, dnrm2

      !------->          ABBREVIATIONS
      !
      !     ne    : number of occupied states
      !     nv    : number of g-components in eigenstate
      !     cv=z  : wavefunction in g-space (reciprocal space)
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

      stateRadius = MAXVAL(ABS(lapw%rk(:,:)))
      stateRadius = stateRadius + SQRT(DOT_PRODUCT(kpts%bk(:,ikpt),kpts%bk(:,ikpt)))
      IF (noco%l_noco) stateRadius = stateRadius + SQRT(DOT_PRODUCT(nococonv%qss(:),nococonv%qss(:)))

      stateFFTRadius = 2.0*stateRadius
      stateFFTExtendedRadius = 2.0*stateRadius

      IF (noco%l_noco.OR.noco%l_soc) THEN
         IF (banddos%dos .OR. banddos%vacdos .OR. input%cdinf.OR.banddos%band) THEN
            stateFFTExtendedRadius = 3.0*stateRadius+0.1
            CALL stepFct%init(cell,sym,stateFFTExtendedRadius+0.001)
            CALL stepFct%putFieldOnGrid(stars, stars%ustep, cell,stateFFTRadius+0.0005)
            CALL stepFct%fillFieldSphereIndexArray(stars, stateFFTRadius+0.0008, fieldSphereIndices)
            CALL fft_interface(3, stepFct%dimensions(:), stepFct%grid, .FALSE., fieldSphereIndices)
         END IF
      END IF

      ALLOCATE (cwk(stars%ng3), ecwk(stars%ng3))

      IF (noco%l_noco) THEN
         CALL state%init(cell,sym,stateFFTExtendedRadius+0.001)
         CALL stateB%init(cell,sym,stateFFTExtendedRadius+0.001)
         DO i = 1, 4
            CALL rhomatGrid(i)%init(cell,sym,stateFFTExtendedRadius+0.001)
            rhomatGrid(i)%grid(:) = CMPLX(0.0,0.0)
         END DO
         CALL rhomatGrid(1)%fillFieldSphereIndexArray(stars, stateFFTRadius+0.0008, fieldSphereIndices)
         IF (noco%l_ss) THEN
            ALLOCATE(stateIndices(lapw%nv(1)))
            ALLOCATE(stateBIndices(lapw%nv(2)))
            CALL state%fillStateIndexArray(lapw,1,stateIndices)
            CALL stateB%fillStateIndexArray(lapw,2,stateBIndices)
         ELSE
            ALLOCATE(stateIndices(lapw%nv(jspin)))
            ALLOCATE(stateBIndices(lapw%nv(jspin)))
            CALL state%fillStateIndexArray(lapw,jspin,stateIndices)
            CALL stateB%fillStateIndexArray(lapw,jspin,stateBIndices)
         ENDIF
      ELSE
         CALL state%init(cell,sym,stateFFTExtendedRadius+0.001)
         CALL chargeDen%init(cell,sym,stateFFTExtendedRadius+0.001)
         chargeDen%grid(:) = CMPLX(0.0,0.0)
         CALL chargeDen%fillFieldSphereIndexArray(stars, stateFFTRadius+0.0008, fieldSphereIndices)
         IF (input%l_f) THEN
            CALL stateDeriv%init(cell,sym,stateFFTExtendedRadius+0.001)
            CALL ekinGrid%init(cell,sym,stateFFTExtendedRadius+0.001)
            ekinGrid%grid(:) = CMPLX(0.0,0.0)
         END IF
         ALLOCATE(stateIndices(lapw%nv(jspin)))
         CALL state%fillStateIndexArray(lapw,jspin,stateIndices)
      ENDIF

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

      IF ((noco%l_noco).AND.(ikpt.LE.fmpi%isize)) THEN
         dos%qis = 0.0
      END IF

      wtf(:ne) = we(:ne)/cell%omtil

      ! LOOP OVER OCCUPIED STATES
      DO nu = 1, ne

         ! Do inverse FFT of state and add related charge density to overall charge density on real-space mesh.
         IF (noco%l_noco) THEN
            forw = .FALSE.
            IF (noco%l_ss) THEN
               CALL state%putComplexStateOnGrid(lapw, 1, zMat%data_c(1:lapw%nv(1),nu))
               CALL stateB%putComplexStateOnGrid(lapw, 2, zMat%data_c(lapw%nv(1) + atoms%nlotot+1:lapw%nv(1) + atoms%nlotot+lapw%nv(2),nu))
            ELSE
               CALL state%putComplexStateOnGrid(lapw, jspin, zMat%data_c(1:lapw%nv(jspin),nu))
               CALL stateB%putComplexStateOnGrid(lapw, jspin, zMat%data_c(lapw%nv(1) + atoms%nlotot+1:lapw%nv(1) + atoms%nlotot+lapw%nv(jspin),nu))
            ENDIF

            CALL fft_interface(3, state%dimensions(:), state%grid, forw, stateIndices)
            CALL fft_interface(3, stateB%dimensions(:), stateB%grid, forw, stateBIndices)

            ! In the non-collinear case rho becomes a hermitian 2x2
            ! matrix (rhomatGrid).
            DO ir = 0, rhomatGrid(1)%gridLength - 1
               rhomatGrid(1)%grid(ir) = rhomatGrid(1)%grid(ir) + wtf(nu) * ABS(state%grid(ir))**2
               rhomatGrid(2)%grid(ir) = rhomatGrid(2)%grid(ir) + wtf(nu) * ABS(stateB%grid(ir))**2
               rhomatGrid(3)%grid(ir) = rhomatGrid(3)%grid(ir) + wtf(nu) * (REAL(state%grid(ir))*REAL(stateB%grid(ir)) + AIMAG(state%grid(ir))*AIMAG(stateB%grid(ir)))
               rhomatGrid(4)%grid(ir) = rhomatGrid(4)%grid(ir) + wtf(nu) * (AIMAG(state%grid(ir))*REAL(stateB%grid(ir)) - REAL(state%grid(ir))*AIMAG(stateB%grid(ir)))
            END DO

            ! In a non-collinear calculation the interstitial charge
            ! cannot be calculated by a simple substraction if the
            ! muffin-tin (and vacuum) charge is know, because the
            ! total charge does not need to be one in each spin-
            ! channel. Thus it has to be calculated explicitly, if
            ! it is needed.
            IF (banddos%dos .OR. banddos%vacdos .OR. input%cdinf.OR.banddos%band) THEN
               DO ir = 0, state%gridLength - 1
                  state%grid(ir) = ABS(state%grid(ir))**2
                  stateB%grid(ir) = ABS(stateB%grid(ir))**2
               END DO

               forw = .TRUE.
               CALL fft_interface(3, state%dimensions(:), state%grid, forw, fieldSphereIndices)
               CALL fft_interface(3, stateB%dimensions(:), stateB%grid, forw, fieldSphereIndices)

               CALL state%takeFieldFromGrid(stars, cwk, stateFFTRadius+0.0005)
               DO istr = 1, stars%ng3
                  cwk(istr) = REAL(stars%nstr(istr))*cwk(istr)
               END DO
               DO istr = 1, stars%ng3_fft
                  CALL pwint(stars, atoms, sym, oneD, cell, istr, x)
                  dos%qis(ev_list(nu), ikpt, 1) = dos%qis(ev_list(nu), ikpt, 1) + REAL(cwk(istr)*x)/cell%omtil
                  dos%qTot(ev_list(nu), ikpt, 1) = dos%qTot(ev_list(nu), ikpt, 1) + REAL(cwk(istr)*x)/cell%omtil
               ENDDO

               CALL stateB%takeFieldFromGrid(stars, cwk, stateFFTRadius+0.0005)
               DO istr = 1, stars%ng3
                  cwk(istr) = REAL(stars%nstr(istr))*cwk(istr)
               END DO
               DO istr = 1, stars%ng3_fft
                  CALL pwint(stars, atoms, sym, oneD, cell, istr, x)
                  dos%qis(ev_list(nu), ikpt, input%jspins) = dos%qis(ev_list(nu), ikpt, input%jspins) + REAL(cwk(istr)*x)/cell%omtil
                  dos%qTot(ev_list(nu), ikpt, input%jspins) = dos%qTot(ev_list(nu), ikpt, input%jspins) + REAL(cwk(istr)*x)/cell%omtil
               ENDDO
            ENDIF

         ELSE

            CALL state%putStateOnGrid(lapw, jSpin, zMat, nu)
            
            forw = .FALSE.
            ! The following FFT is a general complex to complex FFT
            ! For zmat%l_real this should be turned into areal to real FFT at some point
            ! Note: For the moment also no zero-indices for SpFFT provided
            CALL fft_interface(3, state%dimensions(:), state%grid, forw, stateIndices)
            DO ir = 0, chargeDen%gridLength - 1
               chargeDen%grid(ir) = chargeDen%grid(ir) + wtf(nu) * ABS(state%grid(ir))**2
            END DO

            IF (input%l_f) THEN
               DO ir = 0, ekinGrid%gridLength - 1
                  ekinGrid%grid(ir) = ekinGrid%grid(ir) - wtf(nu) * eig(nu) * ABS(state%grid(ir))**2
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
                  CALL fft_interface(3, stateDeriv%dimensions(:), stateDeriv%grid, forw, stateIndices)
                  DO ir = 0, ekinGrid%gridLength - 1
                     ekinGrid%grid(ir) = ekinGrid%grid(ir) + wtf(nu) * 0.5 * ABS(stateDeriv%grid(ir))**2
                  END DO
               END DO
            END IF

            IF (noco%l_soc.AND.input%jspins.EQ.2) THEN
               IF (banddos%dos .OR. banddos%vacdos .OR. input%cdinf.OR.banddos%band) THEN
                  DO ir = 0, state%gridLength - 1
                     state%grid(ir) = conjg(state%grid(ir)) * stepFct%grid(ir) * state%grid(ir)
                  END DO

                  forw = .TRUE.
                  CALL fft_interface(3, state%dimensions(:), state%grid, forw, fieldSphereIndices)

                  cwk = CMPLX(0.0,0.0)
                  CALL state%takeFieldFromGrid(stars, cwk, stateFFTRadius+0.0005)
                  starCharges = CMPLX(0.0,0.0)
                  CALL pwint_all(stars,atoms,sym,oneD,cell,1,stars%ng3,starCharges)
                  starCharges(:) = starCharges(:) * cwk(:) * stars%nstr(:) / cell%omtil
                  dos%qis(ev_list(nu), ikpt, jSpin) = dos%qis(ev_list(nu), ikpt, jSpin) + REAL(SUM(starCharges(:)))
                  dos%qTot(ev_list(nu), ikpt, jSpin) = dos%qTot(ev_list(nu), ikpt, jSpin) + REAL(SUM(starCharges(:)))
               END IF
            END IF
         END IF
      END DO
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
            CALL fft_interface(3, rhomatGrid(idens)%dimensions(:), rhomatGrid(idens)%grid, forw, fieldSphereIndices)
         ELSE
            ! The following FFT is a general complex to complex FFT
            ! For zmat%l_real this should be turned into areal to real FFT at some point
            ! Note: For the moment also no zero-indices for SpFFT provided
            CALL fft_interface(3, chargeDen%dimensions(:), chargeDen%grid, forw, fieldSphereIndices)
            IF (input%l_f) THEN
               CALL fft_interface(3, ekinGrid%dimensions(:), ekinGrid%grid, forw, fieldSphereIndices)
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

         ! add charge density to existing one
         IF (idens .LE. 2) THEN
            ! add to spin-up or -down density (collinear & non-collinear)
            ispin = jspin
            IF (noco%l_noco) ispin = idens
            DO istr = 1, stars%ng3_fft
               den%pw(istr, ispin) = den%pw(istr, ispin) + cwk(istr)
            ENDDO
         ELSE IF (idens .EQ. 3) THEN
            ! add to off-diag. part of density matrix (only non-collinear)
            DO istr = 1, stars%ng3_fft
               den%pw(istr, 3) = den%pw(istr, 3) + cwk(istr)
            ENDDO
         ELSE
            ! add to off-diag. part of density matrix (only non-collinear)
            DO istr = 1, stars%ng3_fft
               den%pw(istr, 3) = den%pw(istr, 3) + CMPLX(0.0, 1.0)*cwk(istr) !TODO: Should be - for no magic minus.
            ENDDO
         ENDIF
      ENDDO

      DEALLOCATE (cwk, ecwk)

      CALL timestop("pwden")

   END SUBROUTINE pwden
END MODULE m_pwden
