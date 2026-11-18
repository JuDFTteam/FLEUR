!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Regression test for the spin rotation in rotMMPmat_dwgn (global/rotMMPmat.f90).
!
! A three-block occupation matrix n(m,m',1:3) stores the physical 2(2l+1) matrix
!    P(up m,up m') = n(m,m',1), P(dn m,dn m') = n(m,m',2),
!    P(dn m,up m') = n(m,m',3), P(up m,dn m') = CONJG(n(m',m,3)).
! rotMMPmat(...,spin_rotation=.TRUE.,real_space_rotation=.FALSE.) intends
!    P -> (su x 1) P (su x 1)^dagger,  su = nococonv%chi(alpha,beta),
! which is computed here independently with full matrix products and compared with the
! routine. Further checks: Hermiticity of the diagonal spin blocks, trace, spectrum,
! forward/inverse round trip, composition of rotations about one axis and the four-block
! branch. INFO lines compare with an emulation of the defective in-place loop (m outer,
! m' inner, read-after-write) and show that inverse=.TRUE. drops alpha from the spin
! rotation (rotMMPmat_angle passes chi(-gamma,-beta)); the latter is not checked here.
!
! Exit status 1 if a check fails; with the defect present, the checks for beta /= 0 fail.
!--------------------------------------------------------------------------------
PROGRAM rotmmpmat_spin_test
   USE m_constants, ONLY: lmaxU_const, ImagUnit
   USE m_rotMMPmat, ONLY: rotMMPmat
   IMPLICIT NONE

   INTEGER, PARAMETER :: lu = lmaxU_const
   REAL,    PARAMETER :: tol = 1.0e-12
   REAL,    PARAMETER :: pi = 3.141592653589793238462643
   INTEGER :: failures, checks, l, ia
   REAL    :: angles(2,6)

   failures = 0
   checks = 0
   CALL init_random()

   ! (alpha,beta): pure z, pure y, the C1 case, the lflip case theta=pi, two general ones
   angles(:,1) = [1.1, 0.0]
   angles(:,2) = [0.0, 0.9]
   angles(:,3) = [0.4, 0.9]
   angles(:,4) = [0.0, pi]
   angles(:,5) = [2.3, 2.0]
   angles(:,6) = [-0.7, 0.3]

   DO l = 2, 3
      DO ia = 1, SIZE(angles,2)
         CALL test_case(l,angles(1,ia),angles(2,ia))
      ENDDO
      CALL test_composition(l,0.5,0.7)
   ENDDO
   CALL test_four_blocks(3)

   WRITE (*,'(a,i0,a,i0)') "rotmmpmat_spin_test: checks = ", checks, ", failures = ", failures
   IF (failures > 0) STOP 1

CONTAINS

   SUBROUTINE test_case(l,alpha,beta)
      INTEGER, INTENT(IN) :: l
      REAL,    INTENT(IN) :: alpha, beta
      COMPLEX, ALLOCATABLE :: v(:,:,:), vrot(:,:,:), vref(:,:,:), vback(:,:,:)
      COMPLEX :: su(2,2)
      CHARACTER(len=40) :: tag
      REAL :: e_in(2*(2*l+1)), e_out(2*(2*l+1))
      REAL :: err, err_lower, err_upper, herm, dtr

      WRITE (tag,'("l=",i1," alpha=",f6.3," beta=",f6.3)') l, alpha, beta
      v = random_occupation(l)
      su = chi(alpha,beta)
      vref = rotate_reference(v,su,l)
      vrot = rotMMPmat(v,alpha,beta,0.0,l,spin_rotation=.TRUE.,real_space_rotation=.FALSE.)

      ! diagnostics of the defective in-place loop (m outer, m' inner); other loop orders differ
      WRITE (*,'(2a,es10.3)') TRIM(tag), ": INFO |rotMMPmat - in-place emulation| = ", &
                             MAXVAL(ABS(vrot-rotate_inplace_emulation(v,su,.FALSE.)))
      IF (ABS(beta) > 1.0e-6) THEN
         WRITE (*,'(2a,es10.3)') TRIM(tag), ": INFO |in-place(m outer) - in-place(m'' outer)| = ", &
                                MAXVAL(ABS(rotate_inplace_emulation(v,su,.FALSE.)-rotate_inplace_emulation(v,su,.TRUE.)))
      ENDIF

      err = MAXVAL(ABS(vrot-vref))
      err_upper = max_part(vrot-vref,l,.FALSE.)
      err_lower = max_part(vrot-vref,l,.TRUE.)
      WRITE (*,'(2a,3es10.3)') TRIM(tag), ": err all / m<=m'' / m>m'' = ", err, err_upper, err_lower
      CALL check(TRIM(tag)//": rotMMPmat = (su x 1) P (su x 1)^dagger", err)

      herm = MAX(MAXVAL(ABS(vrot(:,:,1)-CONJG(TRANSPOSE(vrot(:,:,1))))), &
                 MAXVAL(ABS(vrot(:,:,2)-CONJG(TRANSPOSE(vrot(:,:,2))))))
      CALL check(TRIM(tag)//": spin-diagonal blocks Hermitian", herm)

      dtr = ABS(trace(vrot(:,:,1))+trace(vrot(:,:,2))-trace(v(:,:,1))-trace(v(:,:,2)))
      CALL check(TRIM(tag)//": trace", dtr)

      e_in = hermitian_eigenvalues(physical(v,l))
      e_out = hermitian_eigenvalues(physical(vrot,l))
      CALL check(TRIM(tag)//": spectrum of Hermitian part", MAXVAL(ABS(e_out-e_in)))

      ! inverse=.TRUE. with alpha=0: R^-1 R = 1 (alpha /= 0 is a separate issue, see below)
      IF (ABS(alpha) < 1.0e-12) THEN
         vback = rotMMPmat(vrot,alpha,beta,0.0,l,spin_rotation=.TRUE.,real_space_rotation=.FALSE.,inverse=.TRUE.)
         CALL check(TRIM(tag)//": inverse round trip", MAXVAL(ABS(vback-v)))
      ELSE
         vback = rotMMPmat(vref,alpha,beta,0.0,l,spin_rotation=.TRUE.,real_space_rotation=.FALSE.,inverse=.TRUE.)
         WRITE (*,'(2a,es10.3)') TRIM(tag), ": INFO inverse=.TRUE. applied to the exact rotation, |result - P| = ", &
                                MAXVAL(ABS(vback-v))
         WRITE (*,'(2a,es10.3)') TRIM(tag), ": INFO same with the exact inverse rotate(P, chi(0,-beta)), |diff| = ", &
                                MAXVAL(ABS(rotate_reference(vref,chi(0.0,-beta),l)-vback))
      ENDIF
   END SUBROUTINE test_case

   SUBROUTINE test_composition(l,b1,b2)
      !! rotations about one axis compose: R(0,b2) R(0,b1) = R(0,b1+b2)
      INTEGER, INTENT(IN) :: l
      REAL,    INTENT(IN) :: b1, b2
      COMPLEX, ALLOCATABLE :: v(:,:,:), v1(:,:,:), v12(:,:,:), vdirect(:,:,:)
      CHARACTER(len=40) :: tag

      WRITE (tag,'("l=",i1," composition")') l
      v = random_occupation(l)
      v1 = rotMMPmat(v,0.0,b1,0.0,l,spin_rotation=.TRUE.,real_space_rotation=.FALSE.)
      v12 = rotMMPmat(v1,0.0,b2,0.0,l,spin_rotation=.TRUE.,real_space_rotation=.FALSE.)
      vdirect = rotMMPmat(v,0.0,b1+b2,0.0,l,spin_rotation=.TRUE.,real_space_rotation=.FALSE.)
      CALL check(TRIM(tag)//": R(b2)R(b1) = R(b1+b2)", MAXVAL(ABS(v12-vdirect)))
   END SUBROUTINE test_composition

   SUBROUTINE test_four_blocks(l)
      !! four-block branch: d(1,2) is read from n(m,m',4); with su = 1 nothing may change
      INTEGER, INTENT(IN) :: l
      COMPLEX, ALLOCATABLE :: v3(:,:,:), v4(:,:,:), vrot(:,:,:)
      COMPLEX :: id(2,2)
      INTEGER :: m, mp

      v3 = random_occupation(l)
      ALLOCATE(v4(-lu:lu,-lu:lu,4))
      v4(:,:,1:3) = v3
      DO m = -lu, lu
         DO mp = -lu, lu
            v4(m,mp,4) = CONJG(v3(mp,m,3))
         ENDDO
      ENDDO
      id = CMPLX(0.0,0.0)
      id(1,1) = CMPLX(1.0,0.0)
      id(2,2) = CMPLX(1.0,0.0)
      vrot = rotMMPmat(v4,su=id,spin_rotation=.TRUE.,real_space_rotation=.FALSE.)
      CALL check("four blocks, su = 1: unchanged", MAXVAL(ABS(vrot-v4)))
   END SUBROUTINE test_four_blocks

   FUNCTION random_occupation(l) RESULT(v)
      !! three blocks of a random Hermitian physical matrix, zero outside |m| <= l
      INTEGER, INTENT(IN) :: l
      COMPLEX, ALLOCATABLE :: v(:,:,:)
      COMPLEX :: p(2*(2*l+1),2*(2*l+1))
      REAL    :: re(2*(2*l+1),2*(2*l+1)), im(2*(2*l+1),2*(2*l+1))
      INTEGER :: n

      n = 2*l+1
      CALL RANDOM_NUMBER(re)
      CALL RANDOM_NUMBER(im)
      p = CMPLX(re-0.5,im-0.5)
      p = 0.5*(p+CONJG(TRANSPOSE(p)))
      ALLOCATE(v(-lu:lu,-lu:lu,3),source=CMPLX(0.0,0.0))
      v(-l:l,-l:l,1) = p(1:n,1:n)
      v(-l:l,-l:l,2) = p(n+1:2*n,n+1:2*n)
      v(-l:l,-l:l,3) = p(n+1:2*n,1:n)
   END FUNCTION random_occupation

   FUNCTION physical(v,l) RESULT(p)
      !! physical matrix as FLEUR reads the three blocks; (up,dn) = n(:,:,3)^dagger
      COMPLEX, INTENT(IN) :: v(-lu:,-lu:,:)
      INTEGER, INTENT(IN) :: l
      COMPLEX :: p(2*(2*l+1),2*(2*l+1))
      INTEGER :: n

      n = 2*l+1
      p(1:n,1:n) = v(-l:l,-l:l,1)
      p(n+1:2*n,n+1:2*n) = v(-l:l,-l:l,2)
      p(n+1:2*n,1:n) = v(-l:l,-l:l,3)
      p(1:n,n+1:2*n) = CONJG(TRANSPOSE(v(-l:l,-l:l,3)))
   END FUNCTION physical

   FUNCTION rotate_reference(v,su,l) RESULT(vout)
      !! (su x 1) P (su x 1)^dagger with full matrix products
      COMPLEX, INTENT(IN) :: v(-lu:,-lu:,:)
      COMPLEX, INTENT(IN) :: su(2,2)
      INTEGER, INTENT(IN) :: l
      COMPLEX, ALLOCATABLE :: vout(:,:,:)
      COMPLEX :: p(2*(2*l+1),2*(2*l+1)), u(2*(2*l+1),2*(2*l+1))
      INTEGER :: n, i

      n = 2*l+1
      u = CMPLX(0.0,0.0)
      DO i = 1, n
         u(i,i) = su(1,1)
         u(i,n+i) = su(1,2)
         u(n+i,i) = su(2,1)
         u(n+i,n+i) = su(2,2)
      ENDDO
      p = MATMUL(u,MATMUL(physical(v,l),CONJG(TRANSPOSE(u))))
      ALLOCATE(vout(-lu:lu,-lu:lu,3),source=CMPLX(0.0,0.0))
      vout(-l:l,-l:l,1) = p(1:n,1:n)
      vout(-l:l,-l:l,2) = p(n+1:2*n,n+1:2*n)
      vout(-l:l,-l:l,3) = p(n+1:2*n,1:n)
   END FUNCTION rotate_reference

   FUNCTION rotate_inplace_emulation(v,su,mp_outer) RESULT(out)
      !! the spin loop of rotMMPmat_dwgn, reading d(1,2) from the array being written
      COMPLEX, INTENT(IN) :: v(-lu:,-lu:,:)
      COMPLEX, INTENT(IN) :: su(2,2)
      LOGICAL, INTENT(IN) :: mp_outer
      COMPLEX, ALLOCATABLE :: out(:,:,:)
      COMPLEX :: d(2,2)
      INTEGER :: i, j, m, mp

      out = v
      DO i = -lu, lu
         DO j = -lu, lu
            m = MERGE(j,i,mp_outer)
            mp = MERGE(i,j,mp_outer)
            d(1,1) = out(m,mp,1)
            d(2,2) = out(m,mp,2)
            d(2,1) = out(m,mp,3)
            d(1,2) = CONJG(out(mp,m,3))
            d = MATMUL(su,MATMUL(d,CONJG(TRANSPOSE(su))))
            out(m,mp,1) = d(1,1)
            out(m,mp,2) = d(2,2)
            out(m,mp,3) = d(2,1)
         ENDDO
      ENDDO
   END FUNCTION rotate_inplace_emulation

   FUNCTION chi(alpha,beta) RESULT(su)
      !! nococonv%chi(alpha,beta) (types_nococonv.F90), written out
      REAL, INTENT(IN) :: alpha, beta
      COMPLEX :: su(2,2)

      su(1,1) =  EXP( ImagUnit*alpha/2)*COS(beta/2)
      su(2,1) = -EXP( ImagUnit*alpha/2)*SIN(beta/2)
      su(1,2) =  EXP(-ImagUnit*alpha/2)*SIN(beta/2)
      su(2,2) =  EXP(-ImagUnit*alpha/2)*COS(beta/2)
      su = TRANSPOSE(CONJG(su))
   END FUNCTION chi

   REAL FUNCTION max_part(a,l,lower)
      !! max |a| over m > m' (lower) or m <= m', all blocks
      COMPLEX, INTENT(IN) :: a(-lu:,-lu:,:)
      INTEGER, INTENT(IN) :: l
      LOGICAL, INTENT(IN) :: lower
      INTEGER :: m, mp

      max_part = 0.0
      DO m = -l, l
         DO mp = -l, l
            IF ((m > mp) .EQV. lower) max_part = MAX(max_part,MAXVAL(ABS(a(m,mp,:))))
         ENDDO
      ENDDO
   END FUNCTION max_part

   COMPLEX FUNCTION trace(a)
      COMPLEX, INTENT(IN) :: a(-lu:,-lu:)
      INTEGER :: m

      trace = CMPLX(0.0,0.0)
      DO m = -lu, lu
         trace = trace + a(m,m)
      ENDDO
   END FUNCTION trace

   FUNCTION hermitian_eigenvalues(p) RESULT(e)
      COMPLEX, INTENT(IN) :: p(:,:)
      REAL    :: e(SIZE(p,1))
      COMPLEX :: h(SIZE(p,1),SIZE(p,1)), work(4*SIZE(p,1))
      REAL    :: rwork(3*SIZE(p,1))
      INTEGER :: info

      h = 0.5*(p+CONJG(TRANSPOSE(p)))
      CALL zheev('N','U',SIZE(p,1),h,SIZE(p,1),e,work,SIZE(work),rwork,info)
      IF (info /= 0) e = HUGE(1.0)
   END FUNCTION hermitian_eigenvalues

   SUBROUTINE check(label,err)
      CHARACTER(len=*), INTENT(IN) :: label
      REAL,             INTENT(IN) :: err

      checks = checks + 1
      IF (err > tol) THEN
         failures = failures + 1
         WRITE (*,'(a,a,es10.3)') "FAIL ", label//": ", err
      ENDIF
   END SUBROUTINE check

   SUBROUTINE init_random()
      INTEGER :: n
      INTEGER, ALLOCATABLE :: seed(:)

      CALL RANDOM_SEED(size=n)
      ALLOCATE(seed(n))
      seed = 20261118
      CALL RANDOM_SEED(put=seed)
   END SUBROUTINE init_random

END PROGRAM rotmmpmat_spin_test
