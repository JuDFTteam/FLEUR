!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The single place that knows the on-disk layout of the real-space Wannier
!>  operators O(R) consumed by external transport post-processing. The gauge
!>  rotation and the Fourier transform are operator-independent (melem_ft); what
!>  differs per operator is only the format, selected here by `fmt`:
!>    'hr'      Wannier90 seedname_hr.dat : header + ndegen block + H(R) in eV
!>    'r'       Wannier90 seedname_r.dat  : header + A(R) in Angstrom, 3 components
!>    'bmn'     Wannier90-like            : header + B(R)=<0n|H r|Rm> in eV*Ang, 3 comps
!>    'soc'     R1 R2 R3  i j jj ii  Re Im   (2x2 spinor blocks -> rssocmat.1)
!>    'generic' R1 R2 R3  i j comp   Re Im   (spin -> rspauli.1, orbital -> anglmomrs.*)
!>    'cart2'   R1 R2 R3  i j        then nine (alpha,beta) components on the line
MODULE m_melem_io
   USE m_juDFT
   USE m_constants, ONLY: hartree_to_ev_const
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_write_realspace

CONTAINS

   !> Write O(R) in the requested layout. Rank 0 only (early return elsewhere, so it
   !> is safe to call from all ranks). ndegen is used by the 'hr' layout only.
   SUBROUTINE melem_write_realspace(o_r, irvec, ndegen, nrpts, nw, ncomp, fmt, fname, irank)
      COMPLEX,          INTENT(IN) :: o_r(:, :, :, :)   ! (nw, nw, nrpts, ncomp)
      INTEGER,          INTENT(IN) :: irvec(:, :)       ! (3, nrpts)
      INTEGER,          INTENT(IN) :: ndegen(:)         ! (nrpts) -- 'hr' only
      INTEGER,          INTENT(IN) :: nrpts, nw, ncomp, irank
      CHARACTER(LEN=*), INTENT(IN) :: fmt, fname

      REAL, PARAMETER :: bohr2ang = 0.5291772109
      INTEGER :: irpt, i, j, kk, ii, jj, c, iu

      IF (irank /= 0) RETURN

      OPEN(newunit=iu, file=TRIM(fname), status='replace')
      SELECT CASE (TRIM(fmt))
      CASE ('hr')
         WRITE(iu,'(a)') ' written by FLEUR wannierlib : H(R) in eV, W90 hr format'
         WRITE(iu,'(i12)') nw
         WRITE(iu,'(i12)') nrpts
         c = 0
         DO irpt = 1, nrpts
            WRITE(iu,'(i5)',advance='no') ndegen(irpt); c = c + 1
            IF (MOD(c,15) == 0) WRITE(iu,'(a)') ''
         END DO
         IF (MOD(c,15) /= 0) WRITE(iu,'(a)') ''
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu,'(5i5,2f12.6)') irvec(:,irpt), i, j, &
               hartree_to_ev_const*REAL(o_r(i,j,irpt,1)), hartree_to_ev_const*AIMAG(o_r(i,j,irpt,1))
         END DO; END DO; END DO
      CASE ('r')
         WRITE(iu,'(a)') ' written by FLEUR wannierlib : A(R)=<0n|r|Rm> in Ang, W90 r format'
         WRITE(iu,'(i12)') nw
         WRITE(iu,'(i12)') nrpts
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu,'(5i5,6f12.6)') irvec(:,irpt), i, j, &
               (bohr2ang*REAL(o_r(i,j,irpt,kk)), bohr2ang*AIMAG(o_r(i,j,irpt,kk)), kk=1,3)
         END DO; END DO; END DO
      CASE ('bmn')
         WRITE(iu,'(a)') ' written by FLEUR wannierlib : B(R)=<0n|H r|Rm> in eV*Ang, 3 components'
         WRITE(iu,'(i12)') nw
         WRITE(iu,'(i12)') nrpts
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu,'(5i5,6f12.6)') irvec(:,irpt), i, j, &
               (hartree_to_ev_const*bohr2ang*REAL(o_r(i,j,irpt,kk)), &
                hartree_to_ev_const*bohr2ang*AIMAG(o_r(i,j,irpt,kk)), kk=1,3)
         END DO; END DO; END DO
      CASE ('soc')
         DO irpt = 1, nrpts; DO i = 1, nw; DO j = 1, nw; DO ii = 1, 2; DO jj = 1, 2
            c = (ii - 1)*2 + jj
            WRITE(iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
               irvec(1,irpt), irvec(2,irpt), irvec(3,irpt), i, j, jj, ii, &
               REAL(o_r(i,j,irpt,c)), AIMAG(o_r(i,j,irpt,c))
         END DO; END DO; END DO; END DO; END DO
      CASE ('generic')
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw; DO kk = 1, ncomp
            WRITE(iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
               irvec(1,irpt), irvec(2,irpt), irvec(3,irpt), i, j, kk, &
               REAL(o_r(i,j,irpt,kk)), AIMAG(o_r(i,j,irpt,kk))
         END DO; END DO; END DO; END DO
      CASE ('cart2')
         !> Two Cartesian indices on one line, as the modern-theory post-processing expects:
         !> R, band pair, then the nine (alpha,beta) components in row-major order. ncomp must
         !> be nine; a single-index operator uses 'generic' instead.
         IF (ncomp /= 9) CALL juDFT_error( &
            'melem_write_realspace: the cart2 format carries nine Cartesian components', &
            hint='use the generic format for an operator with one Cartesian index', &
            calledby='melem_write_realspace')
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,18(1x,f20.8))') &
               irvec(1,irpt), irvec(2,irpt), irvec(3,irpt), i, j, &
               (REAL(o_r(i,j,irpt,kk)), AIMAG(o_r(i,j,irpt,kk)), kk = 1, 9)
         END DO; END DO; END DO
      CASE DEFAULT
         CALL juDFT_error('melem_write_realspace: unknown format "'//TRIM(fmt)//'"', &
                          calledby='melem_write_realspace')
      END SELECT
      CLOSE(iu)
   END SUBROUTINE melem_write_realspace

END MODULE m_melem_io
