!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mixvector
   !TODO!!!
   ! LDA+U
   ! Noco (third spin)
#ifdef CPP_MPI
   use mpi
#endif
   USE m_types
   IMPLICIT NONE

   PRIVATE
   !Here we store the pointers used for metric
   TYPE(t_stars), POINTER  :: stars
   TYPE(t_cell), POINTER   :: cell
   TYPE(t_sphhar), POINTER :: sphhar
   TYPE(t_atoms), POINTER  :: atoms => NULL()
   TYPE(t_sym), POINTER    :: sym => NULL()
   INTEGER                :: jspins, nvac
   LOGICAL                :: l_noco, invs, invs2, l_mtnocopot, l_spinoffd_ldau
   INTEGER                :: pw_length !The shape of the local arrays
   INTEGER                :: pw_start(3) = 0, pw_stop(3) !First and last index for spin
   INTEGER                :: mt_length, mt_length_g
   INTEGER                :: mt_start(3) = 0, mt_stop(3) !First and last index for spin
   INTEGER                :: vac_length, vac_length_g
   INTEGER                :: vac_start(3) = 0, vac_stop(3) !First and last index for spin
   INTEGER                :: misc_length = 0, misc_length_g
   INTEGER                :: misc_start(3) = 0, misc_stop(3) !First and last index for spin
   INTEGER                :: mix_mpi_comm !Communicator for all PEs doing mixing
   LOGICAL                :: spin_here(3) = .TRUE.
   LOGICAL                :: pw_here = .TRUE.
   LOGICAL                :: mt_here = .TRUE.
   LOGICAL                :: vac_here = .TRUE.
   LOGICAL                :: misc_here = .TRUE.
   INTEGER                :: mt_rank = 0
   INTEGER                :: mt_size = 1
   LOGICAL                :: l_pot = .FALSE. !Is this a potential?
   REAL, ALLOCATABLE       :: g_mt(:), g_vac(:), g_misc(:)

   TYPE, PUBLIC:: t_mixvector
      REAL, ALLOCATABLE       :: vec_pw(:)
      REAL, ALLOCATABLE       :: vec_mt(:)
      REAL, ALLOCATABLE       :: vec_vac(:)
      REAL, ALLOCATABLE       :: vec_misc(:)
   CONTAINS
      PROCEDURE :: alloc => mixvector_alloc
      PROCEDURE :: from_density => mixvector_from_density
      PROCEDURE :: to_density => mixvector_to_density
      PROCEDURE :: apply_metric => mixvector_metric
      PROCEDURE :: multiply_dot_mask
      PROCEDURE :: dfpt_multiply_dot_mask
      PROCEDURE :: read_unformatted
      PROCEDURE :: write_unformatted
      PROCEDURE :: allocated => mixvector_allocated
   END TYPE t_mixvector

   INTERFACE OPERATOR(*)
      MODULE PROCEDURE multiply_scalar
      MODULE PROCEDURE multiply_scalar_spin
   END INTERFACE OPERATOR(*)
   INTERFACE OPERATOR(+)
      MODULE PROCEDURE add_vectors
   END INTERFACE OPERATOR(+)
   INTERFACE OPERATOR(-)
      MODULE PROCEDURE subtract_vectors
   END INTERFACE OPERATOR(-)
   INTERFACE OPERATOR(.dot.)
      MODULE PROCEDURE multiply_dot
   END INTERFACE OPERATOR(.dot.)

   PUBLIC :: OPERATOR(+), OPERATOR(-), OPERATOR(*), OPERATOR(.dot.)
   PUBLIC :: mixvector_init, mixvector_reset

CONTAINS

   SUBROUTINE READ_unformatted(this, unit)
      IMPLICIT NONE
      CLASS(t_mixvector), INTENT(INOUT)::this
      INTEGER, INTENT(IN)::unit
      call timestart("read_mixing")
      CALL this%alloc()
      IF (pw_here) READ (unit) this%vec_pw
      IF (mt_here) READ (unit) this%vec_mt
      IF (vac_here) READ (unit) this%vec_vac
      IF (misc_here) READ (unit) this%vec_misc
      call timestop("read_mixing")
   END SUBROUTINE READ_unformatted

   SUBROUTINE write_unformatted(this, unit)
      IMPLICIT NONE
      CLASS(t_mixvector), INTENT(IN)::this
      INTEGER, INTENT(IN)::unit
      call timestart("write_mixing")
      IF (pw_here) WRITE (unit) this%vec_pw
      IF (mt_here) WRITE (unit) this%vec_mt
      IF (vac_here) WRITE (unit) this%vec_vac
      IF (misc_here) WRITE (unit) this%vec_misc
      call timestop("write_mixing")
   END SUBROUTINE write_unformatted

   SUBROUTINE mixvector_reset()
      IMPLICIT NONE
      atoms => NULL()
      sym => NULL()
      IF (ALLOCATED(g_mt)) DEALLOCATE (g_mt)
      IF (ALLOCATED(g_vac)) DEALLOCATE (g_vac)
      IF (ALLOCATED(g_misc)) DEALLOCATE (g_misc)
      !restore defaults
      pw_start = 0
      mt_start = 0
      vac_start = 0
      misc_length = 0
      misc_start = 0
      spin_here = .TRUE.
      pw_here = .TRUE.
      mt_here = .TRUE.
      vac_here = .TRUE.
      misc_here = .TRUE.
      mt_rank = 0
      mt_size = 1
      l_pot = .FALSE. !Is this a potential?
   END SUBROUTINE mixvector_reset

   SUBROUTINE mixvector_from_density(vec, den, swapspin, denIm)
      USE m_types
      IMPLICIT NONE
      CLASS(t_mixvector), INTENT(INOUT)    :: vec
      TYPE(t_potden), INTENT(inout)    :: Den
      LOGICAL, INTENT(IN), OPTIONAL         :: swapspin
      TYPE(t_potden), INTENT(INOUT), OPTIONAL :: denIm
      INTEGER:: js, ii, n, l, iv, j, mmpSize
      CALL den%DISTRIBUTE(mix_mpi_comm)
      IF (PRESENT(denIm)) CALL denIm%DISTRIBUTE(mix_mpi_comm)
      DO js = 1, MERGE(jspins, 3,.NOT. l_noco)
         j = js
         IF (PRESENT(swapspin)) THEN
            IF (swapspin .AND. js < 3) j = MERGE(1, 2, js == 2)
         ENDIF
         IF (spin_here(js)) THEN
            !PW part
            IF (pw_here) THEN
               vec%vec_pw(pw_start(js):pw_start(js) + stars%ng3 - 1) = REAL(den%pw(:, j))
               IF ((.NOT. sym%invs) .OR. (js == 3).OR.PRESENT(denIm)) THEN
                  vec%vec_pw(pw_start(js) + stars%ng3:pw_start(js) + 2*stars%ng3 - 1) = AIMAG(den%pw(:, j))
               ENDIF
               IF ((js == 3).AND.PRESENT(denIm)) THEN
                  vec%vec_pw(pw_start(js) + 2*stars%ng3:pw_start(js) + 3*stars%ng3 - 1) =  REAL(den%pw(:, 4))
                  vec%vec_pw(pw_start(js) + 3*stars%ng3:pw_start(js) + 4*stars%ng3 - 1) = AIMAG(den%pw(:, 4))
               END IF
            ENDIF
            IF (vac_here) THEN
               !This PE stores vac-data
               ii = vac_start(js) - 1
               DO iv = 1, nvac
                  vec%vec_vac(ii + 1:ii + SIZE(den%vacz, 1)) = den%vacz(:, iv, j)
                  ii = ii + SIZE(den%vacz, 1)
                  vec%vec_vac(ii + 1:ii + SIZE(den%vacxy(:, :, iv, js))) = RESHAPE(REAL(den%vacxy(:, :, iv, j)), (/SIZE(den%vacxy(:, :, iv, j))/))
                  ii = ii + SIZE(den%vacxy(:, :, iv, j))
                  IF ((.NOT. sym%invs2) .OR. (js == 3)) THEN
                     vec%vec_vac(ii + 1:ii + SIZE(den%vacxy(:, :, iv, j))) = RESHAPE(AIMAG(den%vacxy(:, :, iv, j)), (/SIZE(den%vacxy(:, :, iv, j))/))
                     ii = ii + SIZE(den%vacxy(:, :, iv, j))
                  ENDIF
                  IF (js > 2) THEN
                     vec%vec_vac(ii + 1:ii + SIZE(den%vacz, 1)) = den%vacz(:, iv, 4)
                     ii = ii + SIZE(den%vacz, 1)
                  ENDIF
               ENDDO
            ENDIF
            IF (mt_here .AND. (js < 3 .OR. l_mtnocopot)) THEN
               !This PE stores some(or all) MT data
               ii = mt_start(js) - 1
               IF (.NOT.PRESENT(denIm)) THEN
                  DO n = mt_rank + 1, atoms%ntype, mt_size
                     DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                        vec%vec_mt(ii + 1:ii + atoms%jri(n)) = den%mt(:atoms%jri(n), l, n, j)
                        ii = ii + atoms%jri(n)
                     ENDDO
                  ENDDO
                  IF (js == 3) THEN !Imaginary part
                     DO n = mt_rank + 1, atoms%ntype, mt_size
                        DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                           vec%vec_mt(ii + 1:ii + atoms%jri(n)) = den%mt(:atoms%jri(n), l, n, 4)
                           ii = ii + atoms%jri(n)
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE ! DFPT mixing
                  DO n = mt_rank + 1, atoms%ntype, mt_size
                     DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                        vec%vec_mt(ii + 1:ii + atoms%jri(n)) = den%mt(:atoms%jri(n), l, n, j)
                        ii = ii + atoms%jri(n)
                     END DO
                  END DO
                  DO n = mt_rank + 1, atoms%ntype, mt_size
                     DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                        vec%vec_mt(ii + 1:ii + atoms%jri(n)) = denIm%mt(:atoms%jri(n), l, n, j)
                        ii = ii + atoms%jri(n)
                     END DO
                  END DO
                  IF (js == 3) THEN !Imaginary part
                     DO n = mt_rank + 1, atoms%ntype, mt_size
                        DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                           vec%vec_mt(ii + 1:ii + atoms%jri(n)) = den%mt(:atoms%jri(n), l, n, 4)
                           ii = ii + atoms%jri(n)
                        END DO
                     END DO
                     DO n = mt_rank + 1, atoms%ntype, mt_size
                        DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                           vec%vec_mt(ii + 1:ii + atoms%jri(n)) = denIm%mt(:atoms%jri(n), l, n, 4)
                           ii = ii + atoms%jri(n)
                        END DO
                     END DO
                  END IF
               END IF
            ENDIF
            IF (misc_here .AND. (js < 3 .OR. l_spinoffd_ldau)) THEN
               mmpSize = SIZE(den%mmpMat(:, :, 1:atoms%n_u, j))
               vec%vec_misc(misc_start(js):misc_start(js) + mmpSize - 1) = RESHAPE(REAL(den%mmpMat(:, :, 1:atoms%n_u, j)), (/mmpSize/))
               vec%vec_misc(misc_start(js) + mmpSize:misc_start(js) + 2*mmpSize - 1) = RESHAPE(AIMAG(den%mmpMat(:, :, 1:atoms%n_u, j)), (/mmpSize/))
            END IF
         END IF
      END DO

   END SUBROUTINE mixvector_from_density

   SUBROUTINE mixvector_to_density(vec, den, denIm)
      USE m_types
      IMPLICIT NONE
      CLASS(t_mixvector), INTENT(IN)    :: vec
      TYPE(t_potden), INTENT(INOUT) :: Den
      TYPE(t_potden), INTENT(INOUT), OPTIONAL :: denIm
      INTEGER:: js, ii, n, l, iv, mmpSize
      LOGICAL :: l_dfpt

      l_dfpt = PRESENT(denIm)

      DO js = 1, MERGE(jspins, 3,.NOT. l_noco)
         IF (spin_here(js)) THEN
            !PW part
            IF (pw_here) THEN
               IF (sym%invs .AND. js < 3 .AND. .NOT. l_dfpt) THEN
                  den%pw(:, js) = vec%vec_pw(pw_start(js):pw_start(js) + stars%ng3 - 1)
               ELSE
                  den%pw(:, js) = CMPLX(vec%vec_pw(pw_start(js):pw_start(js) + stars%ng3 - 1), vec%vec_pw(pw_start(js) + stars%ng3:pw_start(js) + 2*stars%ng3 - 1))
                  IF (l_dfpt.AND.js==3) THEN
                     den%pw(:, 4) = CMPLX(vec%vec_pw(pw_start(js) + 2*stars%ng3:pw_start(js) + 3*stars%ng3 - 1), vec%vec_pw(pw_start(js) + 3*stars%ng3:pw_start(js) + 4*stars%ng3 - 1))
                  END IF
               ENDIF
            ENDIF
            IF (mt_here .AND. (js < 3 .OR. l_mtnocopot)) THEN
               !This PE stores some(or all) MT data
               ii = mt_start(js)
               DO n = mt_rank + 1, atoms%ntype, mt_size
                  DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                     den%mt(:atoms%jri(n), l, n, js) = vec%vec_mt(ii:ii + atoms%jri(n) - 1)
                     ii = ii + atoms%jri(n)
                  ENDDO
               ENDDO
               IF (l_dfpt) THEN
                  DO n = mt_rank + 1, atoms%ntype, mt_size
                     DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                        denIm%mt(:atoms%jri(n), l, n, js) = vec%vec_mt(ii:ii + atoms%jri(n) - 1)
                        ii = ii + atoms%jri(n)
                     ENDDO
                  ENDDO
               END IF
               IF (js == 3) THEN !Imaginary part comes as 4th spin
                  DO n = mt_rank + 1, atoms%ntype, mt_size
                     DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                        den%mt(:atoms%jri(n), l, n, 4) = vec%vec_mt(ii:ii + atoms%jri(n) - 1)
                        ii = ii + atoms%jri(n)
                     ENDDO
                  ENDDO
                  IF (l_dfpt) THEN
                     DO n = mt_rank + 1, atoms%ntype, mt_size
                        DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
                           denIm%mt(:atoms%jri(n), l, n, 4) = vec%vec_mt(ii:ii + atoms%jri(n) - 1)
                           ii = ii + atoms%jri(n)
                        ENDDO
                     ENDDO
                  END IF
               ENDIF
            ENDIF
            IF (vac_here) THEN
               !This PE stores vac-data
               ii = vac_start(js) - 1
               DO iv = 1, nvac
                  den%vacz(:, iv, js) = vec%vec_vac(ii + 1:ii + SIZE(den%vacz, 1))
                  ii = ii + SIZE(den%vacz, 1)
                  IF (sym%invs2 .AND. js < 3) THEN
                     den%vacxy(:, :, iv, js) = RESHAPE(vec%vec_vac(ii + 1:ii + SIZE(den%vacxy(:, :, iv, js))), SHAPE(den%vacxy(:, :, iv, js)))
                     ii = ii + SIZE(den%vacxy(:, :, iv, js))
                  ELSE
                     den%vacxy(:, :, iv, js) = RESHAPE(CMPLX(vec%vec_vac(ii + 1:ii + SIZE(den%vacxy(:, :, iv, js))), &
                                                             vec%vec_vac(ii + SIZE(den%vacxy(:, :, iv, js)) + 1:ii + 2*SIZE(den%vacxy(:, :, iv, js)))), &
                                                       SHAPE(den%vacxy(:, :, iv, js)))
                     ii = ii + 2*SIZE(den%vacxy(:, :, iv, js))
                  ENDIF
                  IF (js > 2) THEN
                     den%vacz(:, iv, 4) = vec%vec_vac(ii + 1:ii + SIZE(den%vacz, 1))
                     ii = ii + SIZE(den%vacz, 1)
                  ENDIF
               ENDDO
            ENDIF
            IF (misc_here .AND. (js < 3 .OR. l_spinoffd_ldau)) THEN
               mmpSize = SIZE(den%mmpMat(:, :, 1:atoms%n_u, js))
               den%mmpMat(:, :, 1:atoms%n_u, js) = RESHAPE(CMPLX(vec%vec_misc(misc_start(js):misc_start(js) + mmpSize - 1), &
                                                                 vec%vec_misc(misc_start(js) + mmpSize:misc_start(js) + 2*mmpSize - 1)), &
                                                           SHAPE(den%mmpMat(:, :, 1:atoms%n_u, js)))
            END IF
         END IF
      ENDDO
      CALL den%collect(mix_mpi_comm)

   END SUBROUTINE mixvector_to_density

   FUNCTION mixvector_metric(vec,l_dfpt) RESULT(mvec)
      USE m_types
      USE m_convol
      IMPLICIT NONE
      CLASS(t_mixvector), INTENT(IN) :: vec
      LOGICAL,            INTENT(IN) :: l_dfpt

      TYPE(t_mixvector)              :: mvec

      INTEGER:: js, ii, n, l, iv
      COMPLEX, ALLOCATABLE::pw(:), pw_w(:)
      call timestart("metric")
      mvec = vec
      IF (pw_here) ALLOCATE (pw(stars%ng3), pw_w(stars%ng3))

      DO js = 1, MERGE(jspins, 3,.NOT. l_noco)
         IF (spin_here(js)) THEN
            !PW part
            IF (pw_here) THEN
               !Put back on g-grid and use convol
               IF (sym%invs .AND. js < 3 .AND. .NOT. l_dfpt) THEN
                  pw(:) = vec%vec_pw(pw_start(js):pw_start(js) + stars%ng3 - 1)
               ELSE
                  pw(:) = CMPLX(vec%vec_pw(pw_start(js):pw_start(js) + stars%ng3 - 1), vec%vec_pw(pw_start(js) + stars%ng3:pw_start(js) + 2*stars%ng3 - 1))
               ENDIF
               CALL convol(stars, pw_w, pw)
               pw_w = pw_w*cell%omtil
               mvec%vec_pw(pw_start(js):pw_start(js) + stars%ng3 - 1) = REAL(pw_w)
               IF ((.NOT. sym%invs) .OR. (js == 3) .OR. l_dfpt) THEN
                  mvec%vec_pw(pw_start(js) + stars%ng3:pw_start(js) + 2*stars%ng3 - 1) = AIMAG(pw_w)
               ENDIF
               IF ((js == 3) .AND. l_dfpt) THEN
                  pw(:) = CMPLX(vec%vec_pw(pw_start(js) + 2*stars%ng3:pw_start(js) + 3*stars%ng3 - 1), vec%vec_pw(pw_start(js) + 3*stars%ng3:pw_start(js) + 4*stars%ng3 - 1))
                  CALL convol(stars, pw_w, pw)
                  pw_w = pw_w*cell%omtil
                  mvec%vec_pw(pw_start(js) + 2*stars%ng3:pw_start(js) + 3*stars%ng3 - 1) =  REAL(pw_w)
                  mvec%vec_pw(pw_start(js) + 3*stars%ng3:pw_start(js) + 4*stars%ng3 - 1) = AIMAG(pw_w)
               END IF
            ENDIF
            IF (mt_here .AND. (js < 3 .OR. l_mtnocopot)) THEN
               !This PE stores some(or all) MT data
               IF (.NOT.l_dfpt) THEN
                  mvec%vec_mt(mt_start(js):mt_start(js) + SIZE(g_mt) - 1) = g_mt*vec%vec_mt(mt_start(js):mt_start(js) + SIZE(g_mt) - 1)
                  IF (js == 3) THEN
                     !Here we have a the imaginary part as well
                     mvec%vec_mt(mt_start(js) + SIZE(g_mt):mt_stop(js)) = g_mt*vec%vec_mt(mt_start(js) + SIZE(g_mt):mt_stop(js))
                  ENDIF
               ELSE
                  mvec%vec_mt(mt_start(js):mt_start(js) + SIZE(g_mt) - 1) = g_mt*vec%vec_mt(mt_start(js):mt_start(js) + SIZE(g_mt) - 1)
                  mvec%vec_mt(mt_start(js) + SIZE(g_mt):mt_start(js) + 2*SIZE(g_mt) - 1) = g_mt*vec%vec_mt(mt_start(js) + SIZE(g_mt):mt_start(js) + 2*SIZE(g_mt) - 1)
                  IF (js == 3) THEN
                     mvec%vec_mt(mt_start(js) + 2*SIZE(g_mt):mt_start(js) + 3*SIZE(g_mt) - 1) = g_mt*vec%vec_mt(mt_start(js) + 2*SIZE(g_mt):mt_start(js) + 3*SIZE(g_mt) - 1)
                     mvec%vec_mt(mt_start(js) + 3*SIZE(g_mt):mt_start(js) + 4*SIZE(g_mt) - 1) = g_mt*vec%vec_mt(mt_start(js) + 3*SIZE(g_mt):mt_start(js) + 4*SIZE(g_mt) - 1)
                  ENDIF
               END IF
            ENDIF
            IF (vac_here) THEN
               mvec%vec_vac(vac_start(js):vac_start(js) + SIZE(g_vac) - 1) = g_vac*vec%vec_vac(vac_start(js):vac_start(js) + SIZE(g_vac) - 1)
               IF (js == 3) THEN !We have some extra data that corresponds to first part of metric
                  mvec%vec_vac(vac_start(js) + SIZE(g_vac):vac_stop(js)) = g_vac(:vac_stop(js) - vac_start(js) - SIZE(g_vac) + 1)*vec%vec_vac(vac_start(js) + SIZE(g_vac):vac_stop(js))
               ENDIF
            ENDIF
            IF (misc_here .AND. (js < 3 .OR. l_spinoffd_ldau)) THEN
               mvec%vec_misc(misc_start(js):misc_stop(js)) = g_misc*vec%vec_misc(misc_start(js):misc_stop(js))
            END IF
         ENDIF
      END DO
      call timestop("metric")
   END FUNCTION mixvector_metric

   SUBROUTINE init_metric(vacuum, stars)
      USE m_metrz0
      IMPLICIT NONE
      !
      TYPE(t_vacuum), INTENT(in) :: vacuum
      TYPE(t_stars),  INTENT(in) :: stars

      INTEGER:: i, n, l, j, ivac, iz, iv2c, k2, iv2
      REAL:: dxn, dxn2, dxn4, dvol, volnstr2
      REAL, ALLOCATABLE:: wght(:)

      IF (mt_here) THEN
         !This PE stores some(or all) MT data
         ALLOCATE (g_mt(mt_length_g))
         i = 0
         DO n = mt_rank + 1, atoms%ntype, mt_size
            dxn = atoms%neq(n)*atoms%dx(n)/3.0
            dxn2 = 2.0*dxn
            dxn4 = 4.0*dxn
            DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
               i = i + 1
               g_mt(i) = dxn/atoms%rmsh(1, n)
               IF (.NOT. l_pot) THEN
                  DO j = 2, atoms%jri(n) - 1, 2
                     i = i + 2
                     g_mt(i - 1) = dxn4/atoms%rmsh(j, n)
                     g_mt(i) = dxn2/atoms%rmsh(j + 1, n)
                  END DO
                  ! CHANGE JR 96/12/01
                  ! take care when jri(n) is even
                  i = i + 1 - MOD(atoms%jri(n), 2)
                  g_mt(i) = dxn/atoms%rmsh(atoms%jri(n), n)
               ELSE
                  ! for the potential multiply by r^4
                  DO j = 2, atoms%jri(n) - 1, 2
                     i = i + 2
                     g_mt(i - 1) = dxn4*atoms%rmsh(j, n)**3
                     g_mt(i) = dxn2*atoms%rmsh(j + 1, n)**3
                  END DO
                  i = i + 1 - MOD(atoms%jri(n), 2)
                  g_mt(i) = dxn*atoms%rmsh(atoms%jri(n), n)**3
               END IF
            END DO
         END DO
      ENDIF
      i = 0
      IF (vac_here) THEN
         iv2 = 2
         IF (sym%invs2) iv2 = 1

         ALLOCATE (g_vac(vac_length_g), wght(vacuum%nmzd))
         g_vac(:) = 0.0
         dvol = cell%area*vacuum%delz
         ! nvac=1 if (zrfs.or.invs)
         IF (vacuum%nvac .EQ. 1) dvol = dvol + dvol
         DO ivac = 1, vacuum%nvac
            ! G||=0 components
            !
            ! use 7-point simpson integration in accordance to intgz0.f
            ! calculate weights for integration
            CALL metr_z0(vacuum%nmz, wght)
            DO iz = 1, vacuum%nmz
               i = i + 1
               !
               g_vac(i) = wght(iz)*dvol
               !
            END DO
            ! G||.ne.0 components
            !
            ! calculate weights for integration
            CALL metr_z0(vacuum%nmzxy, wght)
            DO iv2c = 1, iv2
               DO k2 = 1, stars%ng2 - 1
                  !
                  volnstr2 = dvol*stars%nstr2(k2)
                  DO iz = 1, vacuum%nmzxy
                     i = i + 1
                     g_vac(i) = wght(iz)*volnstr2
                  END DO
                  !
               END DO
            END DO
         END DO
      END IF
      IF (misc_here) THEN
         ALLOCATE (g_misc(misc_length_g))
         g_misc = 1.0
      END IF

   END SUBROUTINE init_metric

   SUBROUTINE init_storage_mpi(comm_mpi)
      IMPLICIT NONE
      INTEGER, INTENT(in):: comm_mpi
      INTEGER      :: irank, isize, err, js, new_comm
      mix_mpi_comm = comm_mpi
#ifdef CPP_MPI

      CALL mpi_comm_rank(comm_mpi, irank, err)
      CALL mpi_comm_size(comm_mpi, isize, err)

      IF (isize == 1) RETURN !No parallelization
      js = MERGE(jspins, 3,.NOT. l_noco)!distribute spins
      js = MIN(js, isize)
      CALL judft_comm_split(comm_mpi, MOD(irank, js), irank, new_comm)
      spin_here = (/MOD(irank, js) == 0, MOD(irank, js) == 1, (isize == 2 .AND. irank == 0) .OR. MOD(irank, js) == 2/)

      CALL mpi_comm_rank(new_comm, irank, err)
      CALL mpi_comm_size(new_comm, isize, err)
      CALL mpi_comm_free(new_comm, err)

      !Now distribute data
      IF (isize == 1) RETURN !No further parallelism
      !Split off the pw-part
      pw_here = (irank == 0)
      mt_here = (irank > 0)
      vac_here = vac_here .AND. (irank > 0)
      misc_here = misc_here .AND. (irank > 0)
      isize = isize - 1
      irank = irank - 1
      mt_rank = irank
      mt_size = isize
      IF (isize == 1 .OR. irank < 0) RETURN !No further parallelism
      IF (vac_here .OR. misc_here) THEN !split off-vacuum&misc part
         vac_here = vac_here .AND. (irank == 0)
         misc_here = misc_here .AND. (irank == 0)
         mt_here = (irank > 0)
         isize = isize - 1
         irank = irank - 1
      ENDIF
      mt_rank = irank
      mt_size = isize
#endif
   END SUBROUTINE init_storage_mpi

   SUBROUTINE mixvector_init(comm_mpi, l_densitymatrix,   input, vacuum, noco, stars_i, cell_i, sphhar_i, atoms_i, sym_i, l_dfpt)
      USE m_types
      IMPLICIT NONE
      INTEGER, INTENT(IN)               :: comm_mpi
      LOGICAL, INTENT(IN)               :: l_densitymatrix

      TYPE(t_input), INTENT(IN)         :: input
      TYPE(t_vacuum), INTENT(IN), TARGET :: vacuum
      TYPE(t_noco), INTENT(IN)          :: noco
      TYPE(t_stars), INTENT(IN), TARGET  :: stars_i
      TYPE(t_cell), INTENT(IN), TARGET   :: cell_i
      TYPE(t_sphhar), INTENT(IN), TARGET :: sphhar_i
      TYPE(t_atoms), INTENT(IN), TARGET  :: atoms_i
      TYPE(t_sym), INTENT(IN), TARGET    :: sym_i

      LOGICAL, INTENT(IN) :: l_dfpt

      INTEGER :: js, n, len

      !Store pointers to data-types
      IF (ASSOCIATED(atoms)) RETURN !was done before...
      jspins = input%jspins
      nvac = vacuum%nvac
      l_noco = noco%l_noco
      l_mtnocopot = any(noco%l_unrestrictMT)
      l_spinoffd_ldau = any(noco%l_unrestrictMT).OR.any(noco%l_spinoffd_ldau)
      stars => stars_i; cell => cell_i; sphhar => sphhar_i; atoms => atoms_i; sym => sym_i

      vac_here = input%film
      misc_here = l_densitymatrix
      CALL init_storage_mpi(comm_mpi)

      pw_length = 0; mt_length = 0; vac_length = 0; misc_length = 0
      mt_length_g = 0; vac_length_g = 0; misc_length_g = 0
      DO js = 1, MERGE(jspins, 3,.NOT. l_noco)
         IF (spin_here(js)) THEN
            !Now calculate the length of the vectors
            IF (pw_here) THEN
               pw_start(js) = pw_length + 1
               IF (sym%invs .AND. js < 3 .AND. .NOT. l_dfpt) THEN
                  pw_length = pw_length + stars%ng3
               ELSE
                  pw_length = pw_length + 2*stars%ng3
               ENDIF
               IF (l_dfpt.AND.js==3) pw_length = pw_length + 2*stars%ng3
            ENDIF
            pw_stop(js) = pw_length
            IF (mt_here) THEN
               IF (js < 3 .OR. any(noco%l_unrestrictMT)) mt_start(js) = mt_length + 1
               len = 0
               !This PE stores some(or all) MT data
               DO n = mt_rank + 1, atoms%ntype, mt_size
                  IF (l_dfpt) THEN
                     len = len + 2*(sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1)) + 1)*atoms%jri(n)
                  ELSE
                     len = len + (sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1)) + 1)*atoms%jri(n)
                  END IF
               ENDDO
               mt_length_g = MAX(len, mt_length_g)
               IF (l_dfpt) mt_length_g = mt_length_g / 2
               IF (js == 3) THEN
                  !need to store imaginary part as well...
                  DO n = mt_rank + 1, atoms%ntype, mt_size
                     IF (l_dfpt) THEN
                        len = len + 2*(sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1)) + 1)*atoms%jri(n)
                     ELSE
                        len = len + (sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1)) + 1)*atoms%jri(n)
                     END IF
                  ENDDO
               ENDIF
               IF (js < 3 .OR. any(noco%l_unrestrictMT)) mt_length = mt_length + len
               mt_stop(js) = mt_length
            END IF
            IF (vac_here) THEN
               !This PE stores vac-data
               vac_start(js) = vac_length + 1
               len = 0
               IF (sym%invs2 .AND. js < 3) THEN
                  len = len + vacuum%nmzxyd*(stars%ng2 - 1)*vacuum%nvac + vacuum%nmzd*vacuum%nvac
               ELSE
                  len = len + 2*vacuum%nmzxyd*(stars%ng2 - 1)*vacuum%nvac + vacuum%nmzd*vacuum%nvac
               ENDIF
               vac_length_g = MAX(vac_length_g, len)
               IF (js == 3) len = len + vacuum%nmzd*vacuum%nvac !Offdiagnal potential is complex
               vac_length = vac_length + len
               vac_stop(js) = vac_length
            ENDIF
            IF (misc_here .AND. (js < 3 .OR. l_spinoffd_ldau)) THEN
               len = 7*7*2*atoms%n_u
               misc_start(js) = misc_length + 1
               misc_length = misc_length + len
               misc_stop(js) = misc_length
               misc_length_g = MAX(len, misc_length_g)
            END IF
         END IF
      END DO
      CALL init_metric(vacuum, stars)
   END SUBROUTINE mixvector_init
   SUBROUTINE mixvector_alloc(vec)
      IMPLICIT NONE
      CLASS(t_mixvector), INTENT(OUT)    :: vec
      ALLOCATE (vec%vec_pw(pw_length))
      ALLOCATE (vec%vec_mt(mt_length))
      ALLOCATE (vec%vec_vac(vac_length))
      ALLOCATE (vec%vec_misc(misc_length))
   END SUBROUTINE mixvector_alloc

   FUNCTION multiply_scalar(scalar, vec) RESULT(vecout)
      TYPE(t_mixvector), INTENT(IN)::vec
      REAL, INTENT(IN)             ::scalar
      TYPE(t_mixvector)           ::vecout

      vecout = vec
      vecout%vec_pw = vecout%vec_pw*scalar
      vecout%vec_mt = vecout%vec_mt*scalar
      vecout%vec_vac = vecout%vec_vac*scalar
      vecout%vec_misc = vecout%vec_misc*scalar
   END FUNCTION multiply_scalar

   FUNCTION multiply_scalar_spin(scalar, vec) RESULT(vecout)
      TYPE(t_mixvector), INTENT(IN)::vec
      REAL, INTENT(IN)             ::scalar(:)
      TYPE(t_mixvector)           ::vecout

      INTEGER:: js
      REAL:: fac

      vecout = vec
      DO js = 1, MERGE(jspins, 3,.NOT. l_noco)
         IF (SIZE(scalar) < js) THEN
            fac = 0.0
         ELSE
            fac = scalar(js)
         ENDIF
         IF (pw_start(js) > 0) vecout%vec_pw(pw_start(js):pw_stop(js)) = vecout%vec_pw(pw_start(js):pw_stop(js))*fac
         IF (mt_start(js) > 0) vecout%vec_mt(mt_start(js):mt_stop(js)) = vecout%vec_mt(mt_start(js):mt_stop(js))*fac
         IF (vac_start(js) > 0) vecout%vec_vac(vac_start(js):vac_stop(js)) = vecout%vec_vac(vac_start(js):vac_stop(js))*fac
         IF (misc_start(js) > 0) vecout%vec_misc(misc_start(js):misc_stop(js)) = vecout%vec_misc(misc_start(js):misc_stop(js))*fac
      END DO
   END FUNCTION multiply_scalar_spin

   FUNCTION add_vectors(vec1, vec2) RESULT(vecout)
      TYPE(t_mixvector), INTENT(IN)::vec1, vec2
      TYPE(t_mixvector)           ::vecout

      vecout = vec1
      vecout%vec_pw = vecout%vec_pw + vec2%vec_pw
      vecout%vec_mt = vecout%vec_mt + vec2%vec_mt
      vecout%vec_vac = vecout%vec_vac + vec2%vec_vac
      vecout%vec_misc = vecout%vec_misc + vec2%vec_misc
   END FUNCTION add_vectors

   FUNCTION subtract_vectors(vec1, vec2) RESULT(vecout)
      TYPE(t_mixvector), INTENT(IN)::vec1, vec2
      TYPE(t_mixvector)           ::vecout

      vecout = vec1
      vecout%vec_pw = vecout%vec_pw - vec2%vec_pw
      vecout%vec_mt = vecout%vec_mt - vec2%vec_mt
      vecout%vec_vac = vecout%vec_vac - vec2%vec_vac
      vecout%vec_misc = vecout%vec_misc - vec2%vec_misc
   END FUNCTION subtract_vectors

   FUNCTION multiply_dot(vec1, vec2) RESULT(dprod)
      TYPE(t_mixvector), INTENT(IN)::vec1, vec2
      REAL                        ::dprod, dprod_tmp
      INTEGER                     ::ierr
      dprod = DOT_PRODUCT(vec1%vec_pw, vec2%vec_pw)
      dprod = dprod + DOT_PRODUCT(vec1%vec_mt, vec2%vec_mt)
      dprod = dprod + DOT_PRODUCT(vec1%vec_vac, vec2%vec_vac)
      dprod = dprod + DOT_PRODUCT(vec1%vec_misc, vec2%vec_misc)
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(dprod, dprod_tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mix_mpi_comm, ierr)
      dprod = dprod_tmp
#endif
   END FUNCTION multiply_dot

   FUNCTION multiply_dot_mask(vec1, vec2, mask, spin) RESULT(dprod)
      CLASS(t_mixvector), INTENT(IN)::vec1
      TYPE(t_mixvector), INTENT(IN)::vec2
      LOGICAL, INTENT(IN)          ::mask(4)
      INTEGER, INTENT(IN)          ::spin
      REAL                        ::dprod, dprod_tmp

      INTEGER:: js, ierr

      dprod = 0.0

      DO js = 1, 3
         IF (mask(1) .AND. (spin == js .OR. spin == 0) .AND. pw_start(js) > 0) &
            dprod = dprod + DOT_PRODUCT(vec1%vec_pw(pw_start(js):pw_stop(js)), &
                                        vec2%vec_pw(pw_start(js):pw_stop(js)))
         IF (mask(2) .AND. (spin == js .OR. spin == 0) .AND. mt_start(js) > 0) &
            dprod = dprod + DOT_PRODUCT(vec1%vec_mt(mt_start(js):mt_stop(js)), &
                                        vec2%vec_mt(mt_start(js):mt_stop(js)))
         IF (mask(3) .AND. (spin == js .OR. spin == 0) .AND. vac_start(js) > 0) &
            dprod = dprod + DOT_PRODUCT(vec1%vec_vac(vac_start(js):vac_stop(js)), &
                                        vec2%vec_vac(vac_start(js):vac_stop(js)))
         IF (mask(4) .AND. (spin == js .OR. spin == 0) .AND. misc_start(js) > 0) &
            dprod = dprod + DOT_PRODUCT(vec1%vec_misc(misc_start(js):misc_stop(js)), &
                                        vec2%vec_misc(misc_start(js):misc_stop(js)))
      ENDDO

#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(dprod, dprod_tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mix_mpi_comm, ierr)
      dprod = dprod_tmp
#endif
   END FUNCTION multiply_dot_mask

   SUBROUTINE dfpt_multiply_dot_mask(vec1, vec2, mask, spin, dprod1, dprod2)
      CLASS(t_mixvector), INTENT(IN)::vec1
      TYPE(t_mixvector),  INTENT(IN)::vec2

      LOGICAL, INTENT(IN)    :: mask(2)
      INTEGER, INTENT(IN)    :: spin
      REAL,    INTENT(INOUT) :: dprod1(2)

      REAL, OPTIONAL, INTENT(INOUT) :: dprod2(2)

      REAL :: dprod1_tmp(2), dprod2_tmp(2)
      INTEGER:: js, ierr

      dprod1 = 0.0
      IF (PRESENT(dprod2)) dprod2 = 0.0

      DO js = 1, 2
         IF (mask(1) .AND. (spin == js) .AND. pw_start(js) > 0) THEN
            dprod1(1) = dprod1(1) + DOT_PRODUCT(vec1%vec_pw(pw_start(js):pw_stop(js)/2), &
                                                vec2%vec_pw(pw_start(js):pw_stop(js)/2))
            dprod1(2) = dprod1(2) + DOT_PRODUCT(vec1%vec_pw(pw_stop(js)/2+1:pw_stop(js)), &
                                                vec2%vec_pw(pw_stop(js)/2+1:pw_stop(js)))
         END IF
         IF (mask(2) .AND. (spin == js) .AND. mt_start(js) > 0) THEN
            dprod1(1) = dprod1(1) + DOT_PRODUCT(vec1%vec_mt(mt_start(js):mt_stop(js)/2), &
                                                vec2%vec_mt(mt_start(js):mt_stop(js)/2))
            dprod1(2) = dprod1(2) + DOT_PRODUCT(vec1%vec_mt(mt_stop(js)/2+1:mt_stop(js)), &
                                                vec2%vec_mt(mt_stop(js)/2+1:mt_stop(js)))
         END IF
      END DO

      IF (js==3.AND.PRESENT(dprod2)) THEN
         IF (mask(1) .AND. pw_start(js) > 0) THEN
            dprod1(1) = dprod1(1) + DOT_PRODUCT(vec1%vec_pw(pw_start(js):pw_stop(js)/4), &
                                                vec2%vec_pw(pw_start(js):pw_stop(js)/4))
            dprod1(2) = dprod1(2) + DOT_PRODUCT(vec1%vec_pw(pw_stop(js)/4+1:pw_stop(js)/2), &
                                                vec2%vec_pw(pw_stop(js)/4+1:pw_stop(js)/2))
            dprod2(1) = dprod2(1) + DOT_PRODUCT(vec1%vec_pw(pw_stop(js)/2+1:3*pw_stop(js)/4), &
                                                vec2%vec_pw(pw_stop(js)/2+1:3*pw_stop(js)/4))
            dprod2(2) = dprod2(2) + DOT_PRODUCT(vec1%vec_pw(3*pw_stop(js)/4+1:pw_stop(js)), &
                                                vec2%vec_pw(3*pw_stop(js)/4+1:pw_stop(js)))
         END IF
         IF (mask(2) .AND. pw_start(js) > 0) THEN
            dprod1(1) = dprod1(1) + DOT_PRODUCT(vec1%vec_mt(mt_start(js):mt_stop(js)/4), &
                                                vec2%vec_mt(mt_start(js):mt_stop(js)/4))
            dprod1(2) = dprod1(2) + DOT_PRODUCT(vec1%vec_mt(mt_stop(js)/4+1:mt_stop(js)/2), &
                                                vec2%vec_mt(mt_stop(js)/4+1:mt_stop(js)/2))
            dprod2(1) = dprod2(1) + DOT_PRODUCT(vec1%vec_mt(mt_stop(js)/2+1:3*mt_stop(js)/4), &
                                                vec2%vec_mt(mt_stop(js)/2+1:3*mt_stop(js)/4))
            dprod2(2) = dprod2(2) + DOT_PRODUCT(vec1%vec_mt(3*mt_stop(js)/4+1:mt_stop(js)), &
                                                vec2%vec_mt(3*mt_stop(js)/4+1:mt_stop(js)))
         END IF
      END IF

#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(dprod1, dprod1_tmp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, mix_mpi_comm, ierr)
      dprod1 = dprod1_tmp
      IF (PRESENT(dprod2)) THEN
         CALL MPI_ALLREDUCE(dprod2, dprod2_tmp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, mix_mpi_comm, ierr)
         dprod2 = dprod2_tmp
      END IF
#endif
   END SUBROUTINE dfpt_multiply_dot_mask

   FUNCTION mixvector_allocated(self) RESULT(l_array)
      IMPLICIT NONE
      CLASS(t_mixvector), INTENT(in) :: self
      LOGICAL, ALLOCATABLE :: l_array(:)

      l_array = [ALLOCATED(self%vec_pw), &
                 ALLOCATED(self%vec_mt), &
                 ALLOCATED(self%vec_vac), &
                 ALLOCATED(self%vec_misc)]
   END FUNCTION mixvector_allocated
END MODULE m_types_mixvector
