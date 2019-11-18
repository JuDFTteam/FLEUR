!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_trafo

CONTAINS

   SUBROUTINE waveftrafo_symm(cmt_out, z_out, cmt, l_real, z_r, z_c, bandi, ndb, &
                              nk, iop, atoms, mpbasis, hybrid, kpts, sym, &
                              jsp, dimension, cell, lapw)

      USE m_constants
      USE m_wrapper
      USE m_types
      USE m_juDFT
      IMPLICIT NONE

      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_mpbasis), INTENT(IN)     :: mpbasis
      TYPE(t_hybrid), INTENT(IN)      :: hybrid
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_cell), INTENT(IN)        :: cell
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_lapw), INTENT(IN)        :: lapw

!     - scalars -
      INTEGER, INTENT(IN)      :: nk, jsp, ndb
      INTEGER, INTENT(IN)      ::  bandi, iop

!     - arrays -
      COMPLEX, INTENT(IN)      ::  cmt(:,:,:)
      LOGICAL, INTENT(IN)      ::  l_real
      REAL, INTENT(IN)         ::  z_r(:,:)
      COMPLEX, INTENT(IN)      ::  z_c(:,:)
      COMPLEX, INTENT(OUT)     ::  cmt_out(hybrid%maxlmindx, atoms%nat, ndb)
      COMPLEX, INTENT(OUT)     ::  z_out(lapw%nv(jsp), ndb)

!     - local -

!     - scalars -
      INTEGER                 ::  iatom, iatom1, iiatom, itype, igpt, igpt1, ieq, ieq1, iiop
      INTEGER                 ::  i, l, n, nn, lm0, lm1, lm2, m1, m2
      COMPLEX                 ::  cdum
      COMPLEX, PARAMETER       ::  img = (0.0, 1.0)

!     - arrays -
      REAL                    ::  rrot(3, 3), invrrot(3, 3)
      INTEGER                 ::  g(3), g1(3)
      REAL                    ::  tau1(3), rtaual(3), rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  cmthlp(2*atoms%lmaxd + 1)
      LOGICAL                 ::  trs


      if (l_real) THEN
         rrot    = transpose(1.0 *  sym%mrot(:, :, sym%invtab(iop)))
         invrrot = transpose(1.0 * sym%mrot(:, :, iop))
         trans   = sym%tau(:, iop)
      else
         IF (iop <= sym%nop) THEN
            trs     = .false.
            rrot    = transpose(1.0 * sym%mrot(:, :, sym%invtab(iop)))
            invrrot = transpose(1.0 * sym%mrot(:, :, iop))
            trans   = sym%tau(:, iop)
         ELSE
            trs     = .true.
            iiop    = iop - sym%nop
            rrot    = -transpose(1.0 * sym%mrot(:, :, sym%invtab(iiop)))
            invrrot = -transpose(1.0 * sym%mrot(:, :, iiop))
            trans   = sym%tau(:, iiop)
         END IF
      endif

      rkpt = matmul(rrot, kpts%bkf(:, nk))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g1 = nint(rkpt - rkpthlp)

! MT coefficients
      cmt_out = 0
      iatom = 0
      iiatom = 0

      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            iatom1 = hybrid%map(iatom, iop)
            tau1 = hybrid%tvec(:, iatom, iop)

            cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt, tau1))

            lm0 = 0
            DO l = 0, atoms%lmax(itype)
               nn = mpbasis%num_radfun_per_l(l, itype)
               DO n = 1, nn
                  lm1 = lm0 + n
                  lm2 = lm0 + n + 2*l*nn
                  DO i = 1, ndb
                     if (l_real) THEN
                        cmt_out(lm1:lm2:nn, iatom1, i) = cdum* &
     &                       matmul(cmt(bandi + i - 1, lm1:lm2:nn, iatom),&
     &                       sym%d_wgn(-l:l, -l:l, l, iop))
                     else
                        IF (trs) THEN
                           cmthlp(:2*l + 1) = CONJG(cmt(bandi + i - 1, lm1:lm2:nn, iatom))
                        ELSE
                           cmthlp(:2*l + 1) = cmt(bandi + i - 1, lm1:lm2:nn, iatom)
                        ENDIF
                        cmt_out(lm1:lm2:nn, iatom1, i) = cdum*matmul(cmthlp(:2*l + 1), sym%d_wgn(-l:l, -l:l, l, iop))
                     endif
                  END DO
               END DO
               lm0 = lm2
            END DO
         END DO
         iiatom = iiatom + atoms%neq(itype)
      END DO

! PW coefficients
      z_out = 0

      DO igpt = 1, lapw%nv(jsp)
         g = INT(matmul(invrrot, lapw%gvec(:, igpt, jsp) + g1))
!determine number of g
         igpt1 = 0
         DO i = 1, lapw%nv(jsp)
            IF (maxval(abs(g - lapw%gvec(:, i, jsp))) <= 1E-06) THEN
               igpt1 = i
               EXIT
            END IF
         END DO
         IF (igpt1 == 0) THEN
            call judft_error('wavetrafo_symm: rotated G vector not found')
         END IF
         cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt + lapw%gvec(:, igpt, jsp), trans(:)))
         if (l_real) THEN
            z_out(igpt, 1:ndb) = cdum*z_r(igpt1, bandi:bandi + ndb - 1)
         else
            IF (trs) THEN
               z_out(igpt, 1:ndb) = cdum*CONJG(z_c(igpt1, bandi:bandi + ndb - 1))
            ELSE
               z_out(igpt, 1:ndb) = cdum*z_c(igpt1, bandi:bandi + ndb - 1)
            END IF
         endif
      END DO

   END SUBROUTINE waveftrafo_symm

   SUBROUTINE waveftrafo_genwavf( &
      cmt_out, z_rout, z_cout, cmt, l_real, z_r, z_c, nk, iop, atoms, &
      mpbasis, hybrid, kpts, sym, jsp, nbasfcn, dimension, nbands, &
      cell, lapw_nk, lapw_rkpt, phase)

      use m_juDFT
      USE m_constants
      USE m_wrapper
      USE m_types
      IMPLICIT NONE

      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_mpbasis), INTENT(IN)    :: mpbasis
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_lapw), INTENT(IN)    :: lapw_nk, lapw_rkpt
!     - scalars -
      INTEGER, INTENT(IN)      :: nk, jsp, nbasfcn, nbands
      INTEGER, INTENT(IN)      ::  iop
      LOGICAL                 ::  phase
!     - arrays -
      COMPLEX, INTENT(IN)      ::  cmt(:,:,:)
      LOGICAL, INTENT(IN)      :: l_real
      REAL, INTENT(IN)         ::  z_r(:,:)
      REAL, INTENT(INOUT)      ::  z_rout(:,:)
      COMPLEX, INTENT(IN)      ::  z_c(:,:)
      COMPLEX, INTENT(INOUT)   ::  z_cout(:,:)

      COMPLEX, INTENT(INOUT)  ::  cmt_out(:,:,:)
!        - local -

!     - scalars -
      INTEGER                 ::  itype, iatom, iatom1, iiatom, igpt, igpt1, ieq, ieq1, iiop
      INTEGER                 ::  i, l, n, nn, lm0, lm1, lm2, m1, m2
      COMPLEX                 ::  cdum
      LOGICAL                 ::  trs

!     - arrays -
      INTEGER                 ::  rrot(3, 3), invrrot(3, 3)
      INTEGER                 ::  g(3), g1(3)
      REAL                    ::  tau1(3), rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  zhlp(nbasfcn, dimension%neigd2)
      COMPLEX                 ::  cmthlp(2*atoms%lmaxd + 1)

      if (l_real) THEN
         rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
         invrrot = transpose(sym%mrot(:, :, iop))
         trans = sym%tau(:, iop)
      else
         IF (iop <= sym%nop) THEN
            trs = .false.
            rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
            invrrot = transpose(sym%mrot(:, :, iop))
            trans = sym%tau(:, iop)
         ELSE
! in the case of SOC (l_soc=.true.)
! time reversal symmetry is not valid anymore;
! nsym should thus equal nop
            trs = .true.
            iiop = iop - sym%nop
            rrot = -transpose(sym%mrot(:, :, sym%invtab(iiop)))
            invrrot = -transpose(sym%mrot(:, :, iiop))
            trans = sym%tau(:, iiop)
         END IF
      endif

      rkpt = matmul(rrot, kpts%bkf(:, nk))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g1 = nint(rkpt - rkpthlp)

      ! MT coefficients
      cmt_out = 0
      iatom = 0
      iiatom = 0

      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            iatom1 = hybrid%map(iatom, iop)
            tau1 = hybrid%tvec(:, iatom, iop)

            cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt, tau1))

            lm0 = 0
            DO l = 0, atoms%lmax(itype)
               nn = mpbasis%num_radfun_per_l(l, itype)
               DO n = 1, nn
                  lm1 = lm0 + n
                  lm2 = lm0 + n + 2*l*nn

                  DO i = 1, nbands
                     if (l_real) THEN
                        cmt_out(i, lm1:lm2:nn, iatom1) = cdum*matmul(cmt(i, lm1:lm2:nn, iatom),&
             &                             hybrid%d_wgn2(-l:l, -l:l, l, iop))
                     else
                        IF (trs) THEN
                           cmthlp(:2*l + 1) = conjg(cmt(i, lm1:lm2:nn, iatom))
                        ELSE
                           cmthlp(:2*l + 1) = cmt(i, lm1:lm2:nn, iatom)
                        END IF
                        cmt_out(i, lm1:lm2:nn, iatom1) = cdum*matmul(cmthlp(:2*l + 1), hybrid%d_wgn2(-l:l, -l:l, l, iop))
                     endif

                  END DO
               END DO
               lm0 = lm2
            END DO
         END DO
         iiatom = iiatom + atoms%neq(itype)
      END DO

      ! PW coefficients

      zhlp = 0
      DO igpt = 1, lapw_rkpt%nv(jsp)
         g = matmul(invrrot, lapw_rkpt%gvec(:,igpt,jsp) + g1)
         !determine number of g
         igpt1 = 0
         DO i = 1, lapw_nk%nv(jsp)
            IF (all(abs(g - lapw_nk%gvec(:,i, jsp) ) <= 1E-06)) THEN
               igpt1 = i
               EXIT
            END IF
         END DO
         IF (igpt1 == 0) CYCLE
         cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt + lapw_rkpt%gvec(:,igpt,jsp), trans))
         if (l_real) THEN
            zhlp(igpt, :nbands) = cdum*z_r(igpt1, :nbands)
         else
            IF (trs) THEN
               zhlp(igpt, :nbands) = cdum*conjg(z_c(igpt1, :nbands))
            ELSE
               zhlp(igpt, :nbands) = cdum*z_c(igpt1, :nbands)
            END IF
         endif
      END DO

      ! If phase and inversion-sym. is true,
      ! define the phase such that z_out is real.

      IF (phase) THEN
         DO i = 1, nbands
            if (l_real) THEN
               CALL commonphase(cdum, zhlp(:, i), nbasfcn)

               IF (any(abs(aimag(zhlp(:, i)/cdum)) > 1e-8)) THEN
                  WRITE (*, *) maxval(abs(aimag(zhlp(:, i)/cdum)))
                  WRITE (*, *) zhlp
                  call judft_error('waveftrafo1: Residual imaginary part.')
               END IF
               z_rout(:, i) = real(zhlp(:, i)/cdum)
               cmt_out(i, :, :) = cmt_out(i, :, :)/cdum
            else
               z_cout(:, i) = zhlp(:, i)
            endif
         END DO
      ELSE
         if (l_real) THEN
            z_rout = real(zhlp)
         else
            z_cout = zhlp
         endif
      END IF

   END SUBROUTINE waveftrafo_genwavf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Symmetrizes MT part of input matrix according to inversion symmetry.
   ! This is achieved by a transformation to
   !        1/sqrt(2) * ( exp(ikR) Y_lm(r-R) + (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) )
   ! and                                                                                 if R /=0 or m<0
   !        i/sqrt(2) * ( exp(ikR) Y_lm(r-R) - (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) ) .
   !
   !  or
   !        i*Y_l,0(r)                                                                   if R=0,m=0 and l odd
   ! These functions have the property f(-r)=f(r)* which makes the output matrix real symmetric.
   ! (Array mat is overwritten! )

   SUBROUTINE symmetrize(mat, dim1, dim2, imode, lreal, &
                         atoms, lcutm, maxlcutm, nindxm, sym)
      USE m_types
      IMPLICIT NONE
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_sym), INTENT(IN)     :: sym

!     - scalars -
      INTEGER, INTENT(IN)    ::  imode, dim1, dim2
      INTEGER, INTENT(IN)    :: maxlcutm
      LOGICAL, INTENT(IN)    ::  lreal

!     - arrays -
      INTEGER, INTENT(IN)    :: lcutm(:)
      INTEGER, INTENT(IN)    ::  nindxm(0:maxlcutm, atoms%ntype)
      COMPLEX, INTENT(INOUT) ::  mat(dim1,dim2)

!     -local scalars -
      INTEGER               ::  i, j, itype, ieq, ic, ic1, i1, i2, l, m, n, nn, ifac, ishift
      REAL                  ::  rfac, rdum, rmax
      COMPLEX               ::  img = (0.0, 1.0)

!     - local arrays -
      COMPLEX               ::  carr(max(dim1, dim2)), cfac

      rfac = sqrt(0.5)
      cfac = sqrt(0.5)*img
      ic = 0
      i = 0

      DO itype = 1, atoms%ntype
         nn = sum([((2*l + 1)*nindxm(l, itype), l=0, lcutm(itype))])
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            IF (atoms%invsat(ic) == 0) THEN
! if the structure is inversion-symmetric, but the equivalent atom belongs to a different unit cell
! invsat(atom) = 0, invsatnr(atom) = 0
! but we need invsatnr(atom) = natom
               ic1 = ic
            ELSE
               ic1 = sym%invsatnr(ic)
            END IF
!ic1 = invsatnr(ic)
            IF (ic1 < ic) THEN
               i = i + nn
               CYCLE
            END IF
!     IF( ic1 .lt. ic ) cycle
            DO l = 0, lcutm(itype)
               ifac = -1
               DO m = -l, l
                  ifac = -ifac
                  ishift = (ic1 - ic)*nn - 2*m*nindxm(l, itype)
                  DO n = 1, nindxm(l, itype)
                     i = i + 1
                     j = i + ishift
                     IF (ic1 /= ic .or. m < 0) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           carr(:dim2) = mat(i, :)
                           mat(i, :) = (carr(:dim2) + ifac*mat(j, :))*rfac
                           mat(j, :) = (carr(:dim2) - ifac*mat(j, :))*(-cfac)
                        END IF
                        IF (iand(imode, 2) /= 0) THEN
                           carr(:dim1) = mat(:, i)
                           mat(:, i) = (carr(:dim1) + ifac*mat(:, j))*rfac
                           mat(:, j) = (carr(:dim1) - ifac*mat(:, j))*cfac
                        END IF
                     ELSE IF (m == 0 .and. ifac == -1) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           mat(i, :) = -img*mat(i, :)
                        END IF
                        IF (iand(imode, 2) /= 0) THEN
                           mat(:, i) = img*mat(:, i)
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      IF (lreal) THEN


!     ! Determine common phase factor and devide by it to make the output matrix real.
!     rmax = 0
!     DO i = 1,dim1
!     DO j = 1,dim2
!     rdum = abs(real(mat(i,j)))+abs(aimag(mat(i,j)))
!     IF(rdum.gt.1e-6) THEN
!     cfac = mat(i,j)/abs(mat(i,j))
!     GO TO 1
!     ELSE IF(rdum.gt.rmax) THEN
!     cfac = mat(i,j)/abs(mat(i,j))
!     rmax = rdum
!     END IF
!     END DO
!     END DO
!     IF(1-abs(cfac)   .gt.1e-8) THEN ; mat = 0 ; RETURN ; END IF
!     1      IF(abs(1-cfac**2).gt.1e-8) mat = mat/cfac
!
!     IF(any(abs(aimag(mat)).gt.1e-8)) THEN
!     WRITE(*,*) maxval(aimag(mat))
!     call judft_error('symmetrize: Residual imaginary part. Symmetrization failed.')

! Determine common phase factor and divide by it to make the output matrix real.
         CALL commonphase(cfac, mat, dim1*dim2)
         mat = mat/cfac
         IF (any(abs(aimag(mat)) > 1e-8)) &
     &STOP 'symmetrize: Residual imaginary part. Symmetrization failed.'
      END IF

   END SUBROUTINE symmetrize

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Undoes symmetrization with routine symmetrize.
   SUBROUTINE desymmetrize(mat, dim1, dim2, imode, &
                           atoms, lcutm, maxlcutm, nindxm, sym)

      USE m_types
      IMPLICIT NONE
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      INTEGER, INTENT(IN)      ::  imode, dim1, dim2
      INTEGER, INTENT(IN)      :: maxlcutm

!     - arrays -
      INTEGER, INTENT(IN)      :: lcutm(:)
      INTEGER, INTENT(IN)      ::  nindxm(0:maxlcutm, atoms%ntype)
      COMPLEX, INTENT(INOUT)   ::  mat(dim1,dim2)

!     - local scalars -
      INTEGER                 ::  ifac, i, j, itype, ieq, ic, ic1, i1, i2, l, m, n, nn, ishift
      REAL                    ::  rfac1, rfac2
      COMPLEX                 ::  img = (0.0, 1.0)
!     - local arrays -
      COMPLEX                 ::  carr(max(dim1, dim2))

      rfac1 = sqrt(0.5)
      ic = 0
      i = 0
      DO itype = 1, atoms%ntype
         nn = sum([((2*l + 1)*nindxm(l, itype), l=0, lcutm(itype))])
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            IF (atoms%invsat(ic) == 0) THEN
               ! if the structure is inversion-symmetric, but the equivalent atom belongs to a different unit cell
               ! invsat(atom) = 0, invsatnr(atom) =0
               ! but we need invsatnr(atom) = natom
               ic1 = ic
            ELSE
               ic1 = sym%invsatnr(ic)
            END IF
            !ic1 = invsatnr(ic)
            !IF( ic1 .lt. ic ) cycle
            IF (ic1 < ic) THEN
               i = i + nn
               CYCLE
            END IF
            DO l = 0, lcutm(itype)
               ifac = -1
               DO m = -l, l
                  ifac = -ifac
                  rfac2 = rfac1*ifac
                  ishift = (ic1 - ic)*nn - 2*m*nindxm(l, itype)
                  DO n = 1, nindxm(l, itype)
                     i = i + 1
                     j = i + ishift
                     IF (ic1 /= ic .or. m < 0) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           carr(:dim2) = mat(i, :)
                           mat(i, :) = (carr(:dim2) + img*mat(j, :))*rfac1
                           mat(j, :) = (carr(:dim2) - img*mat(j, :))*rfac2
                        END IF
                        IF (iand(imode, 2) /= 0) THEN
                           carr(:dim1) = mat(:, i)
                           mat(:, i) = (carr(:dim1) - img*mat(:, j))*rfac1
                           mat(:, j) = (carr(:dim1) + img*mat(:, j))*rfac2
                        END IF
                     ELSE IF (m == 0 .and. ifac == -1) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           mat(i, :) = img*mat(i, :)
                        END IF
                        IF (iand(imode, 2) /= 0) THEN
                           mat(:, i) = -img*mat(:, i)
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE desymmetrize

   ! bra_trafo1 rotates cprod at ikpt0(<=> not irreducible k-point) to cprod at ikpt1 (bkp(ikpt0)), which is the
   ! symmetrie equivalent one
   ! isym maps ikpt0 on ikpt1

   SUBROUTINE bra_trafo2( &
      l_real, vecout_r, vecin_r, vecout_c, vecin_c, &
      dim, nobd, nbands, ikpt0, ikpt1, iop, sym, &
      mpbasis, hybrid, kpts, cell, atoms, &
      phase)

      !  ikpt0  ::  parent of ikpt1
      !  iop maps ikpt0 on ikpt1

      USE m_constants
      USE m_dwigner
      USE m_util
      USE m_types
      IMPLICIT NONE
      type(t_mpbasis), intent(in)  :: mpbasis
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      INTEGER, INTENT(IN)      ::  ikpt0, ikpt1, iop, dim, nobd, nbands

!     - arrays -

      LOGICAL, INTENT(IN)      :: l_real

      REAL, INTENT(IN)         ::  vecin_r(:,:,:)
      REAL, INTENT(OUT)        ::  vecout_r(:,:,:)
      COMPLEX, INTENT(IN)      ::  vecin_c(:,:,:)
      COMPLEX, INTENT(OUT)     ::  vecout_c(:,:,:)
      COMPLEX, INTENT(OUT)     ::  phase(:,:)

!          - local -

!     - scalars -
      INTEGER                 ::  nrkpt, rcent, itype, ieq, ic, l, n, i, j, nn, i1, i2, j1, j2, m1, m2, ok
      INTEGER                 ::  igptm, igptm2, igptp, icent1, iiatom, iiop, inviop
      COMPLEX                 ::  cexp, cdum
      COMPLEX, PARAMETER       ::  img = (0.0, 1.0)
!     - arrays -

      INTEGER                 ::  rrot(3, 3), invrot(3, 3)
      INTEGER                 ::  pnt(maxval(mpbasis%num_radbasfn), 0:maxval(hybrid%lcutm1), atoms%nat)
      INTEGER                 ::  g(3), g1(3)
      REAL                    ::  rkpt(3), rkpthlp(3), rtaual(3), trans(3)
      REAL                    ::  arg
      COMPLEX                 ::  dwgn(-maxval(hybrid%lcutm1):maxval(hybrid%lcutm1),&
     &                                 -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1))
!       COMPLEX                 ::  vecin1(dim,nobd,nbands),vecout1(dim,nobd,nbands)
      COMPLEX, ALLOCATABLE    ::  vecin1(:, :, :), vecout1(:, :, :)

      call timestart("bra trafo")

      allocate(vecin1(dim, nobd, nbands), &
     &           vecout1(dim, nobd, nbands), stat=ok)
      IF (ok /= 0) &
     &             call judft_error('bra_trafo2: error allocating vecin1 or vecout1')
      vecin1 = 0; vecout1 = 0

      IF (maxval(hybrid%lcutm1) > atoms%lmaxd) call judft_error('bra_trafo2: maxlcutm > atoms%lmaxd')   ! very improbable case

!     transform back to unsymmetrized product basis in case of inversion symmetry
      if (l_real) THEN
         vecin1 = vecin_r
         DO i = 1, nbands
            DO j = 1, nobd
               CALL desymmetrize(vecin1(:hybrid%nbasp, j, i), hybrid%nbasp, 1, 1, &
                                 atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), mpbasis%num_radbasfn, sym)
            END DO
         END DO
      else
         vecin1 = vecin_c
      endif

      IF (iop <= sym%nop) THEN
         inviop = sym%invtab(iop)
         rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
         invrot = sym%mrot(:, :, sym%invtab(iop))
         trans = sym%tau(:, iop)

         dwgn(-maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1)) &
            = hybrid%d_wgn2(-maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1), inviop)

      ELSE
         iiop = iop - sym%nop
         inviop = sym%invtab(iiop) + sym%nop
         rrot = -transpose(sym%mrot(:, :, sym%invtab(iiop)))
         invrot = sym%mrot(:, :, sym%invtab(iiop))
         trans = sym%tau(:, iiop)

         dwgn(-maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1)) &
            = conjg(hybrid%d_wgn2(-maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1), inviop))

      END IF

      rkpt = matmul(rrot, kpts%bkf(:, ikpt0))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g = nint(rkpthlp - rkpt)

#ifdef CPP_DEBUG
      !test
      nrkpt = 0
      DO i = 1, kpts%nkptf
         IF (maxval(abs(rkpt - kpts%bkf(:, i))) <= 1E-06) THEN
            nrkpt = i
            EXIT
         END IF
      END DO
      IF (nrkpt /= ikpt1) THEN
         PRINT *, ikpt0, ikpt1
         PRINT *, kpts%bkf(:, ikpt1)
         PRINT *, kpts%bkf(:, ikpt0)
         PRINT *, rkpt

         call judft_error('bra_trafo2: rotation failed')
      ENDIF
#endif

!     Define pointer to first mixed-basis functions (with m = -l)
      i = 0
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            DO l = 0, hybrid%lcutm1(itype)
               DO n = 1, mpbasis%num_radbasfn(l, itype)
                  i = i + 1
                  pnt(n, l, ic) = i
               END DO
               i = i + mpbasis%num_radbasfn(l, itype)*2*l
            END DO
         END DO
      END DO

!     Multiplication
      ! MT
      cexp = exp(img*tpi_const*dot_product(kpts%bkf(:, ikpt1) + g, trans(:)))
      ic = 0
      iiatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1

            rcent = hybrid%map(ic, iop)

            cdum = cexp*exp(-img*tpi_const*dot_product(g, atoms%taual(:, rcent)))

            DO l = 0, hybrid%lcutm1(itype)
               nn = mpbasis%num_radbasfn(l, itype)
               DO n = 1, nn

                  i1 = pnt(n, l, ic)
                  i2 = i1 + nn*2*l
                  j1 = pnt(n, l, rcent)
                  j2 = j1 + nn*2*l

                  DO i = 1, nbands
                     DO j = 1, nobd
                        vecout1(i1:i2:nn, j, i) = cdum*matmul(vecin1(j1:j2:nn, j, i), dwgn(-l:l, -l:l, l))
                     END DO
                  END DO

               END DO
            END DO
         END DO
         iiatom = iiatom + atoms%neq(itype)
      END DO

      ! PW
      DO igptm = 1, mpbasis%n_g(ikpt0)
         igptp = mpbasis%gptm_ptr(igptm, ikpt0)
         g1 = matmul(rrot, mpbasis%g(:, igptp)) + g
         igptm2 = 0
         DO i = 1, mpbasis%n_g(ikpt1)
            IF (maxval(abs(g1 - mpbasis%g(:, mpbasis%gptm_ptr(i, ikpt1)))) <= 1E-06) THEN
               igptm2 = i
               EXIT
            END IF
         END DO
         IF (igptm2 == 0) THEN
            WRITE (*, *) ikpt0, ikpt1, g1
            WRITE (*, *) mpbasis%n_g(ikpt0), mpbasis%n_g(ikpt1)
            WRITE (*, *)
            WRITE (*, *) igptp, mpbasis%g(:, igptp)
            WRITE (*, *) g
            WRITE (*, *) rrot
            WRITE (*, *) "Failed tests:", g1
            DO i = 1, mpbasis%n_g(ikpt1)
               WRITE (*, *) mpbasis%g(:, mpbasis%gptm_ptr(i, ikpt1))
            ENDDO
            call judft_error('bra_trafo2: G-point not found in G-point set.')
         END IF
         cdum = exp(img*tpi_const*dot_product(kpts%bkf(:, ikpt1) + g1, trans(:)))

         vecout1(hybrid%nbasp + igptm, :, :) = cdum*vecin1(hybrid%nbasp + igptm2, :, :)
      END DO

      deallocate(vecin1)

      if (l_real) THEN
         DO i = 1, nbands
            DO j = 1, nobd

               CALL symmetrize(vecout1(:, j, i), dim, 1, 1, .false., &
                               atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), mpbasis%num_radbasfn, sym)

               CALL commonphase(phase(j, i), vecout1(:, j, i), dim)
               vecout1(:, j, i) = vecout1(:, j, i)/phase(j, i)
               IF (any(abs(aimag(vecout1(:, j, i))) > 1e-8)) THEN
                  WRITE (*, *) vecout1(:, j, i)
                  call judft_error('bra_trafo2: Residual imaginary part.')
               END IF

            END DO
         END DO
      else
         phase = (1.0, 0.0)
      endif

      if (l_real) THEN
         vecout_r = real(vecout1)
      else
         vecout_c = vecout1
      endif
      deallocate(vecout1)
      call timestop("bra trafo")
   END SUBROUTINE bra_trafo2

   ! Determines common phase factor (with unit norm)
   SUBROUTINE commonphase(cfac, carr, n)
      USE m_juDFT
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: n
      COMPLEX, INTENT(IN)      :: carr(n)
      COMPLEX, INTENT(OUT)     :: cfac
      REAL                    :: rdum, rmax
      INTEGER                 :: i

!       IF( all( abs(carr) .lt. 1E-12 ) ) THEN
!         cfac = 1
!         RETURN
!       END IF

      cfac = 0
      rmax = 0
      DO i = 1, n
         rdum = abs(carr(i))
         IF (rdum > 1e-6) THEN; cfac = carr(i)/rdum; EXIT
         ELSE IF (rdum > rmax) THEN; cfac = carr(i)/rdum; rmax = rdum
         END IF
      END DO
      IF (abs(cfac) < 1e-10 .and. all(abs(carr) > 1e-10)) THEN
         WRITE (999, *) carr
         call judft_error('commonphase: Could not determine common phase factor. (Wrote carr to fort.999)')
      END IF
   END SUBROUTINE commonphase

   SUBROUTINE bramat_trafo( &
      vecout, igptm_out, vecin, igptm_in, ikpt0, iop, writevec, pointer, sym, &
      rrot, invrrot, mpbasis, hybrid, kpts, maxlcutm, atoms, lcutm, nindxm, maxindxm, dwgn, nbasp, nbasm)

      USE m_constants
      USE m_util
      USE m_types
      IMPLICIT NONE
      type(t_mpbasis), intent(in) :: mpbasis
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars
      INTEGER, INTENT(IN)      ::  ikpt0, igptm_in, iop, maxindxm
      INTEGER, INTENT(IN)      ::  maxlcutm
      INTEGER, INTENT(IN)      ::  nbasp
      LOGICAL, INTENT(IN)      ::  writevec
      INTEGER, INTENT(OUT)     ::  igptm_out
!     - arrays -
      INTEGER, INTENT(IN)      ::  rrot(:,:), invrrot(:,:)
      INTEGER, INTENT(IN)      :: lcutm(atoms%ntype),&
     &                            nindxm(0:maxlcutm, atoms%ntype)
      INTEGER, INTENT(IN)      :: nbasm(:)
      INTEGER, INTENT(IN)      ::  pointer(&
     &                          minval(mpbasis%g(1, :)) - 1:maxval(mpbasis%g(1, :)) + 1,&
     &                          minval(mpbasis%g(2, :)) - 1:maxval(mpbasis%g(2, :)) + 1,&
     &                          minval(mpbasis%g(3, :)) - 1:maxval(mpbasis%g(3, :)) + 1)

      COMPLEX, INTENT(IN)      ::  vecin(:)
      COMPLEX, INTENT(IN)      ::  dwgn(-maxlcutm:maxlcutm,&
     &                                 -maxlcutm:maxlcutm,&
     &                                         0:maxlcutm)
      COMPLEX, INTENT(OUT)     ::  vecout(nbasm(ikpt0))

!     - private scalars -
      INTEGER                 ::  itype, ieq, ic, l, n, i, nn, i1, i2, j1, j2
      INTEGER                 ::  igptm, igptm2, igptp, isym
      INTEGER                 ::  ikpt1, rcent
      LOGICAL                 ::  trs
      COMPLEX, PARAMETER       ::  img = (0.0, 1.0)
      COMPLEX                 ::  cexp, cdum
!     - private arrays -
      INTEGER                 ::  pnt(maxindxm, 0:maxlcutm, atoms%nat), g(3),&
     &                            g1(3), iarr(mpbasis%n_g(ikpt0))
      REAL                    ::  rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  vecin1(nbasm(ikpt0))
      COMPLEX                 ::  carr(mpbasis%n_g(ikpt0))

      IF (iop <= sym%nop) THEN
         isym = iop
         trs = .false.
         trans = sym%tau(:, isym)
      ELSE
         isym = iop - sym%nop
         trs = .true.
         trans = sym%tau(:, isym)
      END IF

      rkpthlp = matmul(rrot, kpts%bkf(:, ikpt0))
      rkpt = kpts%to_first_bz(rkpthlp)
      g = nint(rkpthlp - rkpt)
      !
      ! determine number of rotated k-point bk(:,ikpt) -> ikpt1
      !
      DO i = 1, kpts%nkpt
         IF (maxval(abs(rkpt - kpts%bkf(:, i))) <= 1E-06) THEN
            ikpt1 = i
            EXIT
         END IF
      END DO

      DO igptm = 1, mpbasis%n_g(ikpt1)
         igptp = mpbasis%gptm_ptr(igptm, ikpt1)
         g1 = matmul(invrrot, mpbasis%g(:, igptp) - g)
         igptm2 = pointer(g1(1), g1(2), g1(3))
         IF (igptm2 == igptm_in) THEN
            igptm_out = igptm
            IF (writevec) THEN
               cdum = exp(img*tpi_const*dot_product(kpts%bkf(:, ikpt1) + mpbasis%g(:, igptp), trans))
               EXIT
            ELSE
               RETURN
            END IF
         END IF
      END DO

!     Transform back to unsymmetrized product basis in case of inversion symmetry.
      vecout = vecin(:nbasm(ikpt0))
      if (sym%invs) CALL desymmetrize(vecout, nbasp, 1, 1, &
                                      atoms, lcutm, maxlcutm, nindxm, sym)

!     Right-multiplication
      ! PW
      IF (trs) THEN; vecin1 = cdum*conjg(vecout)
      ELSE; vecin1 = cdum*vecout
      END IF

!     Define pointer to first mixed-basis functions (with m = -l)
      i = 0
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            DO l = 0, lcutm(itype)
               DO n = 1, nindxm(l, itype)
                  i = i + 1
                  pnt(n, l, ic) = i
               END DO
               i = i + nindxm(l, itype)*2*l
            END DO
         END DO
      END DO

!     Left-multiplication
      ! MT
      cexp = exp(-img*tpi_const*dot_product(kpts%bkf(:, ikpt1) + g, trans))
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            rcent = hybrid%map(ic, sym%invtab(isym))
            cdum = cexp*exp(img*tpi_const*dot_product(g, atoms%taual(:, ic))) !rcent)))
            cdum = conjg(cdum)
            DO l = 0, lcutm(itype)
               nn = nindxm(l, itype)
               DO n = 1, nn

                  i1 = pnt(n, l, ic)
                  i2 = i1 + nn*2*l
                  j1 = pnt(n, l, rcent)
                  j2 = j1 + nn*2*l

                  vecout(i1:i2:nn) = cdum*matmul(dwgn(-l:l, -l:l, l), vecin1(j1:j2:nn))

               END DO
            END DO
         END DO
      END DO

      ! PW
      DO igptm = 1, mpbasis%n_g(ikpt1)
         igptp = mpbasis%gptm_ptr(igptm, ikpt1)
         g1 = matmul(invrrot, mpbasis%g(:, igptp) - g)
         iarr(igptm) = pointer(g1(1), g1(2), g1(3))
         carr(igptm) = exp(-img*tpi_const*dot_product(kpts%bkf(:, ikpt1) + mpbasis%g(:, igptp), trans))
      END DO
      DO i1 = 1, mpbasis%n_g(ikpt1)
         vecout(nbasp + i1) = carr(i1)*vecin1(nbasp + iarr(i1))
      END DO

      ! If inversion symmetry is applicable, symmetrize to make the values real.
      if (sym%invs) CALL symmetrize(vecout, nbasp, 1, 1, .false., &
                                    atoms, lcutm, maxlcutm, nindxm, sym)

   END SUBROUTINE bramat_trafo

END MODULE m_trafo
