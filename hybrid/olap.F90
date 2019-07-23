MODULE m_olap

CONTAINS

!     Calculates plane-wave overlap matrix olap defined by GPT(1:3,1:NGPT).
!     (Muffin-tin spheres are cut out.)
!     olap_pw calculates full overlap matrix

   SUBROUTINE olap_pw(olap, gpt, ngpt, atoms, cell)

      USE m_constants
      USE m_types
      IMPLICIT NONE
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      INTEGER, INTENT(IN)       :: ngpt
!     - arrays -
      INTEGER, INTENT(IN)       :: gpt(3, ngpt)
      TYPE(t_mat)              :: olap
!     - local -
      INTEGER                  :: i, j, itype, icent, ineq
      REAL                     :: g, r, fgr
      COMPLEX, PARAMETER        :: img = (0.0, 1.0)
      INTEGER                  :: dg(3)

      DO i = 1, ngpt
         DO j = 1, i
            dg = gpt(:, j) - gpt(:, i)
            g = gptnorm(dg, cell%bmat)
            IF (g == 0) THEN
               DO itype = 1, atoms%ntype
                  r = atoms%rmt(itype)
                  if (olap%l_real) THEN
                     olap%data_r(i, j) = olap%data_r(i, j) - atoms%neq(itype)*fpi_const*r**3/3/cell%omtil
                  else
                     olap%data_c(i, j) = olap%data_c(i, j) - atoms%neq(itype)*fpi_const*r**3/3/cell%omtil
                  endif
               END DO
            ELSE
               icent = 0
               DO itype = 1, atoms%ntype
                  r = g*atoms%rmt(itype)
                  fgr = fpi_const*(sin(r) - r*cos(r))/g**3/cell%omtil
                  DO ineq = 1, atoms%neq(itype)
                     icent = icent + 1
                     if (olap%l_real) THEN
                        olap%data_r(i, j) = olap%data_r(i, j) - fgr*exp(img*tpi_const*dot_product(dg, atoms%taual(:, icent)))
                     else
                        olap%data_c(i, j) = olap%data_c(i, j) - fgr*exp(img*tpi_const*dot_product(dg, atoms%taual(:, icent)))
                     endif
                  END DO
               END DO
            END IF
            if (olap%l_real) THEN
               IF (i == j) olap%data_r(i, j) = olap%data_r(i, j) + 1
               olap%data_r(j, i) = olap%data_r(i, j)
            else
               IF (i == j) olap%data_c(i, j) = olap%data_c(i, j) + 1
               olap%data_c(j, i) = conjg(olap%data_c(i, j))
            endif
         END DO
      END DO

   END SUBROUTINE olap_pw

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     olap_pwp  calculates upper triangular part of overlap matrix

   SUBROUTINE olap_pwp(l_real, olap_r, olap_c, gpt, ngpt, atoms, cell)

      USE m_constants
      USE m_types
      IMPLICIT NONE
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      INTEGER, INTENT(IN)       :: ngpt
!     - arrays -
      INTEGER, INTENT(IN)       :: gpt(3, ngpt)

      LOGICAL, INTENT(IN)       :: l_real
      REAL, INTENT(OUT)         ::  olap_r(ngpt*(ngpt + 1)/2)
      COMPLEX, INTENT(OUT)      ::  olap_c(ngpt*(ngpt + 1)/2)
!     - local -
      INTEGER                  :: i, j, k, itype, icent, ineq
      REAL                     :: g, r, fgr
      COMPLEX, PARAMETER        :: img = (0.0, 1.0)
      INTEGER                  :: dg(3)

      if (l_real) THEN
         k = 0
         DO i = 1, ngpt
            DO j = 1, i
               k = k + 1
               dg = gpt(:, i) - gpt(:, j)
               g = gptnorm(dg, cell%bmat)
               olap_r(k) = 0
               IF (g == 0) THEN
                  DO itype = 1, atoms%ntype
                     r = atoms%rmt(itype)
                     olap_r(k) = olap_r(k) - atoms%neq(itype)*fpi_const*r**3/3/cell%omtil
                  END DO
               ELSE
                  icent = 0
                  DO itype = 1, atoms%ntype
                     r = g*atoms%rmt(itype)
                     fgr = fpi_const*(sin(r) - r*cos(r))/g**3/cell%omtil
                     DO ineq = 1, atoms%neq(itype)
                        icent = icent + 1
                        olap_r(k) = olap_r(k) - fgr* &
             &               exp(img*tpi_const*dot_product(dg, atoms%taual(:, icent)))
                     END DO
                  END DO
               END IF
               IF (i == j) olap_r(k) = olap_r(k) + 1
            END DO
         END DO
      else
         k = 0
         DO i = 1, ngpt
            DO j = 1, i
               k = k + 1
               dg = gpt(:, i) - gpt(:, j)
               g = gptnorm(dg, cell%bmat)
               olap_c(k) = 0
               IF (g == 0) THEN
                  DO itype = 1, atoms%ntype
                     r = atoms%rmt(itype)
                     olap_c(k) = olap_c(k) - atoms%neq(itype)*fpi_const*r**3/3/cell%omtil
                  END DO
               ELSE
                  icent = 0
                  DO itype = 1, atoms%ntype
                     r = g*atoms%rmt(itype)
                     fgr = fpi_const*(sin(r) - r*cos(r))/g**3/cell%omtil
                     DO ineq = 1, atoms%neq(itype)
                        icent = icent + 1
                        olap_c(k) = olap_c(k) - fgr* &
             &               exp(img*tpi_const*dot_product(dg, atoms%taual(:, icent)))
                     END DO
                  END DO
               END IF
               IF (i == j) olap_c(k) = olap_c(k) + 1
            END DO
         END DO

      endif
   END SUBROUTINE olap_pwp

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SUBROUTINE wfolap_init(olappw, olapmt, gpt,&
  &                       atoms, hybrid, cell,&
  &                       bas1, bas2)

      USE m_util, ONLY: intgrf, intgrf_init
      USE m_types
      IMPLICIT NONE
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - arrays -
      INTEGER, INTENT(IN)       :: gpt(:, :)!(3,ngpt)
      REAL, INTENT(IN)         ::  bas1(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype),&
     &                            bas2(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype)
      REAL, INTENT(OUT)         :: olapmt(hybrid%maxindx, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype)
      TYPE(t_mat), INTENT(INOUT):: olappw

!     - local -
      INTEGER                  :: itype, l, nn, n1, n2

      REAL, ALLOCATABLE         :: gridf(:, :)

      CALL intgrf_init(atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, gridf)
      olapmt = 0
      DO itype = 1, atoms%ntype
         DO l = 0, atoms%lmax(itype)
            nn = hybrid%nindx(l, itype)
            DO n2 = 1, nn
               DO n1 = 1, nn!n2
                  !IF( n1 .gt. 2 .or. n2 .gt. 2) CYCLE
                  olapmt(n1, n2, l, itype) = intgrf( &
         &                            bas1(:, n1, l, itype)*bas1(:, n2, l, itype)&
         &                           + bas2(:, n1, l, itype)*bas2(:, n2, l, itype),&
         &                            atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
!               olapmt(n2,n1,l,itype) = olapmt(n1,n2,l,itype)
               END DO
            END DO
         END DO
      END DO

      CALL olap_pw(olappw, gpt, size(gpt, 2), atoms, cell)

   END SUBROUTINE wfolap_init

   FUNCTION wfolap_inv(cmt1, cpw1, cmt2, cpw2, ngpt1, ngpt2, olappw, olapmt, atoms, hybrid)

      USE m_wrapper
      USE m_types
      IMPLICIT NONE
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      COMPLEX                :: wfolap_inv
      INTEGER, INTENT(IN)     :: ngpt1, ngpt2
!     - arrays -
      COMPLEX, INTENT(IN)     :: cmt1(hybrid%maxlmindx, atoms%nat),&
     &                          cmt2(hybrid%maxlmindx, atoms%nat)
      REAL, INTENT(IN)        :: cpw1(ngpt1)
      COMPLEX, INTENT(IN)     :: cpw2(ngpt2)
      REAL, INTENT(IN)        :: olappw(ngpt1, ngpt2)
      REAL, INTENT(IN)        :: olapmt(hybrid%maxindx, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype)
!     - local -
      INTEGER                :: itype, ieq, iatom, l, m, lm, nn
      COMPLEX                :: carr(ngpt1), cdum
      REAL                   :: rarr1(ngpt1), rarr2(ngpt1), rdum1, rdum2

      wfolap_inv = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            lm = 0
            DO l = 0, atoms%lmax(itype)
               DO M = -l, l
                  nn = hybrid%nindx(l, itype)
                  wfolap_inv = wfolap_inv + &
         &                 dot_product(cmt1(lm + 1:lm + nn, iatom),&
         &                               matmul(olapmt(:nn, :nn, l, itype),&
         &                                       cmt2(lm + 1:lm + nn, iatom)))
                  lm = lm + nn
               END DO
            END DO
         END DO
      END DO

      wfolap_inv = wfolap_inv + dot_product(cpw1, matmul(olappw, cpw2))

!       CALL dgemv('N',ngpt1,ngpt2,1.0,olappw,ngpt1,real(cpw2),1,0.0,rarr1,1)
!       CALL dgemv('N',ngpt1,ngpt2,1.0,olappw,ngpt1,aimag(cpw2),1,0.0,rarr2,1)
!
!       rdum1 = dotprod(cpw1,rarr1)
!       rdum2 = dotprod(cpw1,rarr2)
!       cdum  = cmplx( rdum1, rdum2 )

!       wfolap = wfolap + cdum

   END FUNCTION wfolap_inv
   FUNCTION wfolap_noinv(cmt1, cpw1, cmt2, cpw2, ngpt1, ngpt2, olappw, olapmt, atoms, hybrid)

      USE m_wrapper
      USE m_types
      IMPLICIT NONE
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      COMPLEX                :: wfolap_noinv
      INTEGER, INTENT(IN)     :: ngpt1, ngpt2
!     - arrays -
      COMPLEX, INTENT(IN)     :: cmt1(hybrid%maxlmindx, atoms%nat),&
     &                          cmt2(hybrid%maxlmindx, atoms%nat)
      COMPLEX, INTENT(IN)     :: cpw1(ngpt1)
      COMPLEX, INTENT(IN)     :: cpw2(ngpt2)
      COMPLEX, INTENT(IN)     :: olappw(ngpt1, ngpt2)
      REAL, INTENT(IN)        :: olapmt(hybrid%maxindx, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype)
!     - local -
      INTEGER                :: itype, ieq, iatom, l, m, lm, nn
      COMPLEX                :: carr(ngpt1), cdum
      REAL                   :: rarr1(ngpt1), rarr2(ngpt1), rdum1, rdum2

      wfolap_noinv = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            lm = 0
            DO l = 0, atoms%lmax(itype)
               DO M = -l, l
                  nn = hybrid%nindx(l, itype)
                  wfolap_noinv = wfolap_noinv + &
         &                 dot_product(cmt1(lm + 1:lm + nn, iatom),&
         &                               matmul(olapmt(:nn, :nn, l, itype),&
         &                                       cmt2(lm + 1:lm + nn, iatom)))
                  lm = lm + nn
               END DO
            END DO
         END DO
      END DO

      wfolap_noinv = wfolap_noinv + dot_product(cpw1, matmul(olappw, cpw2))

!       CALL dgemv('N',ngpt1,ngpt2,1.0,olappw,ngpt1,real(cpw2),1,0.0,rarr1,1)
!       CALL dgemv('N',ngpt1,ngpt2,1.0,olappw,ngpt1,aimag(cpw2),1,0.0,rarr2,1)
!
!       rdum1 = dotprod(cpw1,rarr1)
!       rdum2 = dotprod(cpw1,rarr2)
!       cdum  = cmplx( rdum1, rdum2 )

!       wfolap = wfolap + cdum

   END FUNCTION wfolap_noinv

   FUNCTION wfolap1(cmt1, cpw1, cmt2, cpw2, ngpt1, ngpt2, olappw, olapmt,&
  &                atoms, hybrid)

      USE m_types
      IMPLICIT NONE

      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     -scalars -
      COMPLEX                :: wfolap1
      INTEGER, INTENT(IN)     :: ngpt1, ngpt2
!     - arrays -
      COMPLEX, INTENT(IN)     :: cmt1(hybrid%maxlmindx, atoms%nat),&
     &                          cmt2(hybrid%maxlmindx, atoms%nat)
#if ( defined(CPP_INVERSION) )
      REAL, INTENT(IN)        :: cpw1(ngpt1), cpw2(ngpt2)
#else
      COMPLEX, INTENT(IN)     :: cpw1(ngpt1), cpw2(ngpt2)
#endif
#if ( defined(CPP_INVERSION) )
      REAL, INTENT(IN)        :: olappw(ngpt1, ngpt2)
#else
      COMPLEX, INTENT(IN)     :: olappw(ngpt1, ngpt2)
#endif
      REAL, INTENT(IN)        :: olapmt(hybrid%maxindx, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype)

!     - local -
      INTEGER                :: itype, ieq, iatom, l, m, lm, nn

      wfolap1 = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            lm = 0
            DO l = 0, atoms%lmax(itype)
               DO M = -l, l
                  nn = hybrid%nindx(l, itype)
                  wfolap1 = wfolap1 + &
         &                  dot_product(cmt1(lm + 1:lm + nn, iatom),&
         &                                matmul(olapmt(:nn, :nn, l, itype),&
         &                                        cmt2(lm + 1:lm + nn, iatom)))
                  lm = lm + nn
               END DO
            END DO
         END DO
      END DO

      wfolap1 = wfolap1 + dot_product(cpw1, matmul(olappw, cpw2))

   END FUNCTION wfolap1

   FUNCTION wfolap2(cmt1, cpw1, cmt2, cpw2, ngpt1, ngpt2, olappw, olapmt,&
  &                atoms, hybrid)
      USE m_types
      IMPLICIT NONE

      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      COMPLEX                :: wfolap2
      INTEGER, INTENT(IN)     :: ngpt1, ngpt2
!     - arrays -
      COMPLEX, INTENT(IN)     :: cmt1(hybrid%maxlmindx, atoms%nat),&
     &                          cmt2(hybrid%maxlmindx, atoms%nat)
! #if ( defined(CPP_INVERSION) )
!       REAL,INTENT(IN)        :: cpw1(ngpt1)
! #else
      COMPLEX, INTENT(IN)     :: cpw1(ngpt1)
! #endif
      COMPLEX, INTENT(IN)     :: cpw2(ngpt2)
#if ( defined(CPP_INVERSION) )
      REAL, INTENT(IN)        :: olappw(ngpt1, ngpt2)
#else
      COMPLEX, INTENT(IN)     :: olappw(ngpt1, ngpt2)
#endif
      REAL, INTENT(IN)        :: olapmt(hybrid%maxindx, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype)
!     - local -
      INTEGER                :: itype, ieq, ic, l, m, lm, nn

      wfolap2 = 0
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            lm = 0
            DO l = 0, atoms%lmax(itype)
               DO M = -l, l
                  nn = hybrid%nindx(l, itype)
                  wfolap2 = wfolap2 + &
         &                 dot_product(cmt1(lm + 1:lm + nn, ic),&
         &                               matmul(olapmt(:nn, :nn, l, itype),&
         &                                       cmt2(lm + 1:lm + nn, ic)))
                  lm = lm + nn
               END DO
            END DO
         END DO
      END DO

      wfolap2 = wfolap2 + dot_product(cpw1, matmul(olappw, cpw2))

   END FUNCTION wfolap2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     Calculates overlap between core and valence wavefunctions

   SUBROUTINE olap_cv(hybrid, kpts, maxlcutc, maxindxc, atoms,&
  &                   lmaxc, lmaxcd, nindxc,&
  &                   core1, core2, bas1, bas2, cmt, dimension,&
  &                   gridf)

      USE m_util, ONLY: intgrf, intgrf_init, chr
      USE m_types
      IMPLICIT NONE
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      INTEGER, INTENT(IN)    :: maxlcutc, maxindxc, lmaxcd

!     - arrays -
      INTEGER, INTENT(IN)    ::  lmaxc(atoms%ntype)
      INTEGER, INTENT(IN)    ::  nindxc(0:maxlcutc, atoms%ntype)
      REAL, INTENT(IN)       ::  core1(atoms%jmtd, 0:lmaxcd, maxindxc, atoms%ntype),&
     &                          core2(atoms%jmtd, 0:lmaxcd, maxindxc, atoms%ntype)
      REAL, INTENT(IN)       ::  bas1(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype),&
     &                          bas2(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype)
      COMPLEX, INTENT(IN)    ::  cmt(dimension%neigd, kpts%nkpt, hybrid%maxlmindx, atoms%nat)

!     - local scalars -
      INTEGER               :: itype, icent, l, m, lm, i, j

!     - local arrays -
      INTEGER, ALLOCATABLE   :: olapcv_loc(:, :, :, :, :)
      REAL, ALLOCATABLE      :: gridf(:, :)
      REAL, ALLOCATABLE      :: olapcb(:)
      REAL, ALLOCATABLE      :: olapcv_avg(:, :, :, :), olapcv_max(:, :, :, :)
      COMPLEX, ALLOCATABLE   :: olapcv(:, :)
      CHARACTER, PARAMETER  :: lchar(0:38) =&
     &          (/'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',&
     &            'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x',&
     &            'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x'/)

!      CALL intgrf_init(ntype,jmtd,jri,dx,rmsh,gridf)
      WRITE (6, '(/A)') 'Overlap <core|basis>'
      ALLOCATE (olapcb(hybrid%maxindx), olapcv(dimension%neigd, kpts%nkpt),&
     &     olapcv_avg(-maxlcutc:maxlcutc, maxindxc, 0:maxlcutc, atoms%ntype),&
     &     olapcv_max(-maxlcutc:maxlcutc, maxindxc, 0:maxlcutc, atoms%ntype),&
     &     olapcv_loc(2, -maxlcutc:maxlcutc, maxindxc, 0:maxlcutc, atoms%ntype))

      DO itype = 1, atoms%ntype
         IF (atoms%ntype > 1) WRITE (6, '(A,I3)') 'Atom type', itype
         DO l = 0, lmaxc(itype)
            IF (l > atoms%lmax(itype)) THEN
               WRITE (*, *) 'l greater then atoms%lmax(itype)'
               EXIT ! very improbable case
            END IF
!          WRITE(6,8001) (lchar(l),i=1,min(3,nindx(l,itype)))
            DO i = 1, nindxc(l, itype)
               WRITE (6, '(I1,A,2X)', advance='no') i + l, lchar(l)
               DO j = 1, hybrid%nindx(l, itype)
                  olapcb(j) = intgrf(core1(:, l, i, itype)*bas1(:, j, l, itype) +&
         &                            core2(:, l, i, itype)*bas2(:, j, l, itype),&
         &                            atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
                  WRITE (6, '(F10.6)', advance='no') olapcb(j)
               ENDDO

               lm = sum((/(hybrid%nindx(j, itype)*(2*j + 1), j=0, l - 1)/))
               icent = sum(atoms%neq(1:itype - 1)) + 1 ! take first of group of equivalent atoms
               DO M = -l, l
                  olapcv = 0
                  DO j = 1, hybrid%nindx(l, itype)
                     lm = lm + 1
                     olapcv(:, :) = olapcv(:, :) + olapcb(j)*cmt(:, :, lm, icent)
                  END DO
                  olapcv_avg(M, i, l, itype) = sqrt(sum(abs(olapcv(:, :))**2)&
         &                                          /kpts%nkpt/dimension%neigd)
                  olapcv_max(M, i, l, itype) = maxval(abs(olapcv(:, :)))
                  olapcv_loc(:, M, i, l, itype) = maxloc(abs(olapcv(:, :)))
               END DO
               WRITE (6, *)

            END DO
         END DO
      END DO

      WRITE (6, '(/A)') 'Average overlap <core|val>'

      DO itype = 1, atoms%ntype
         IF (atoms%ntype > 1) write (6, '(A,I3)') 'Atom type', itype
         DO l = 0, lmaxc(itype)
            DO i = 1, nindxc(l, itype)
               WRITE (6, '(I1,A,2X)', advance='no') i + l, lchar(l)
               WRITE (6, '('//chr(2*l + 1)//'F10.6)') olapcv_avg(-l:l, i, l, itype)
            END DO
         END DO
      END DO

      WRITE (6, '(/A)') 'Maximum overlap <core|val> at (band/kpoint)'
      DO itype = 1, atoms%ntype
         IF (atoms%ntype > 1) write (6, '(A,I3)') 'Atom type', itype
         DO l = 0, lmaxc(itype)
            DO i = 1, nindxc(l, itype)
               WRITE (6, '(I1,A,2X)', advance='no') i + l, lchar(l)
               WRITE (6, '('//chr(2*l + 1)//&
        &               '(F10.6,'' ('',I3.3,''/'',I4.3,'')''))')&
        &                      (olapcv_max(M, i, l, itype),&
        &                      olapcv_loc(:, M, i, l, itype), M=-l, l)
            END DO
         END DO
      END DO

      DEALLOCATE (olapcb, olapcv, olapcv_avg, olapcv_max, olapcv_loc)

   END SUBROUTINE olap_cv

END MODULE m_olap
