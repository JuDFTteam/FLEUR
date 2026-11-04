!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The single place that knows the on-disk layout of the real-space Wannier
!>  operators O(R) consumed by external transport post-processing. The gauge
!>  rotation and the Fourier transform are operator-independent (wgauge_ft); what
!>  differs per operator is only the format, selected here by `fmt`:
!>
!>  Energies are converted and lengths are not, because they reach this point from
!>  different places: the eigenvalues come from FLEUR in Hartree, while the b-mesh the
!>  position operators are built from comes from Wannier90, already in Angstrom.
!>    'hr'      Wannier90 seedname_hr.dat : header + ndegen block + H(R) in eV
!>    'r'       Wannier90 seedname_r.dat  : header + A(R) in Angstrom, 3 components
!>    'bmn'     Wannier90-like            : header + B(R)=<0n|H (r-R)|Rm> in eV*Ang, 3 comps
!>    'spinor2x2' R1 R2 R3  i j jj ii  Re Im (2x2 spinor blocks -> rssocmat.1)
!>    'generic' R1 R2 R3  i j comp   Re Im   (spin -> rspauli.1, orbital -> anglmomrs.*)
!>    'cart2'   R1 R2 R3  i j        then nine (alpha,beta) components on the line
!>    'cart2e'  as cart2, with the energy converted to eV
!>
!>  The first three drop the header, and 'hr' sends its degeneracies to `ndegen`
!>  instead, when the filename is one the transport code reads: that one takes the
!>  matrix with list-directed input from the first line. See wannierlib_export_bare.
MODULE m_wgauge_io
   USE m_juDFT
   USE m_constants, ONLY: hartree_to_ev_const
   USE m_wannierlib_export_paths, ONLY: wannierlib_export_bare, wannierlib_export_file
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wgauge_write_realspace

CONTAINS

   !> Write O(R) in the requested layout. Rank 0 only (early return elsewhere, so it
   !> is safe to call from all ranks). ndegen is used by the 'hr' layout only.
   SUBROUTINE wgauge_write_realspace(o_r, irvec, ndegen, nrpts, nw, ncomp, fmt, fname, irank)
      COMPLEX,          INTENT(IN) :: o_r(:, :, :, :)   ! (nw, nw, nrpts, ncomp)
      INTEGER,          INTENT(IN) :: irvec(:, :)       ! (3, nrpts)
      INTEGER,          INTENT(IN) :: ndegen(:)         ! (nrpts) -- 'hr' only
      INTEGER,          INTENT(IN) :: nrpts, nw, ncomp, irank
      CHARACTER(LEN=*), INTENT(IN) :: fmt, fname

      INTEGER :: irpt, i, j, kk, ii, jj, c, iu, iu2
      LOGICAL :: bare

      IF (irank /= 0) RETURN

      bare = wannierlib_export_bare(fname)

      OPEN(newunit=iu, file=TRIM(fname), status='replace')
      SELECT CASE (TRIM(fmt))
      CASE ('hr')
         IF (bare) THEN
            !> The degeneracies still have to reach the consumer, so they go to the file
            !> it reads them from, in the same fifteen-per-line block they occupy here.
            !> Splitting them off is the point: the width of that block is ceil(nrpts/15),
            !> and a reader that skips a fixed number of header lines instead gets it
            !> right for one nrpts and silently wrong for the next.
            OPEN(newunit=iu2, file=wannierlib_export_file('ndegen'), status='replace')
            c = 0
            DO irpt = 1, nrpts
               WRITE(iu2,'(i5)',advance='no') ndegen(irpt); c = c + 1
               IF (MOD(c,15) == 0) WRITE(iu2,'(a)') ''
            END DO
            !> The terminator goes out unconditionally, which leaves a trailing empty record
            !> when nrpts divides by fifteen. The consumer reads floor(nrpts/15) lines of
            !> fifteen and then one more WITHOUT testing whether anything is left, so on an
            !> exact multiple that last read has a zero-length list and still needs a record
            !> to advance over: without this it meets end-of-file. nrpts is a product of odd
            !> numbers, and an 8x8x1 mesh gives (2*8-1)^2 = 225, so the case is reachable
            !> with the standard wannierisation mesh rather than only in principle.
            WRITE(iu2,'(a)') ''
            CLOSE(iu2)
         ELSE
            WRITE(iu,'(a)') ' written by FLEUR wannierlib : H(R) in eV, W90 hr format'
            WRITE(iu,'(i12)') nw
            WRITE(iu,'(i12)') nrpts
            c = 0
            DO irpt = 1, nrpts
               WRITE(iu,'(i5)',advance='no') ndegen(irpt); c = c + 1
               IF (MOD(c,15) == 0) WRITE(iu,'(a)') ''
            END DO
            IF (MOD(c,15) /= 0) WRITE(iu,'(a)') ''
         END IF
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu,'(5i5,2f12.6)') irvec(:,irpt), i, j, &
               hartree_to_ev_const*REAL(o_r(i,j,irpt,1)), hartree_to_ev_const*AIMAG(o_r(i,j,irpt,1))
         END DO; END DO; END DO
      CASE ('r')
         IF (.NOT. bare) THEN
            WRITE(iu,'(a)') ' written by FLEUR wannierlib : A(R)=<0n|r|Rm> in Ang, W90 r format'
            WRITE(iu,'(i12)') nw
            WRITE(iu,'(i12)') nrpts
         END IF
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu,'(5i5,6f12.6)') irvec(:,irpt), i, j, &
               (REAL(o_r(i,j,irpt,kk)), AIMAG(o_r(i,j,irpt,kk)), kk=1,3)
         END DO; END DO; END DO
      CASE ('bmn')
         IF (.NOT. bare) THEN
            WRITE(iu,'(a)') ' written by FLEUR wannierlib : B(R)=<0n|H (r-R)|Rm> in eV*Ang, 3 components'
            WRITE(iu,'(i12)') nw
            WRITE(iu,'(i12)') nrpts
         END IF
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu,'(5i5,6f12.6)') irvec(:,irpt), i, j, &
               (hartree_to_ev_const*REAL(o_r(i,j,irpt,kk)), &
                hartree_to_ev_const*AIMAG(o_r(i,j,irpt,kk)), kk=1,3)
         END DO; END DO; END DO
      CASE ('spinor2x2')
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
            'wgauge_write_realspace: the cart2 format carries nine Cartesian components', &
            hint='use the generic format for an operator with one Cartesian index', &
            calledby='wgauge_write_realspace')
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,18(1x,f20.8))') &
               irvec(1,irpt), irvec(2,irpt), irvec(3,irpt), i, j, &
               (REAL(o_r(i,j,irpt,kk)), AIMAG(o_r(i,j,irpt,kk)), kk = 1, 9)
         END DO; END DO; END DO
      CASE ('cart2e')
         !> The cart2 layout, for an operator that carries energy as well. The conversion
         !> is the only difference, and it happens here rather than in the caller because
         !> this is where the unit contract lives.
         IF (ncomp /= 9) CALL juDFT_error( &
            'wgauge_write_realspace: the cart2e format carries nine Cartesian components', &
            hint='use the generic format for an operator with one Cartesian index', &
            calledby='wgauge_write_realspace')
         DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw
            WRITE(iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,18(1x,f20.8))') &
               irvec(1,irpt), irvec(2,irpt), irvec(3,irpt), i, j, &
               (hartree_to_ev_const*REAL(o_r(i,j,irpt,kk)), &
                hartree_to_ev_const*AIMAG(o_r(i,j,irpt,kk)), kk = 1, 9)
         END DO; END DO; END DO
      CASE DEFAULT
         CALL juDFT_error('wgauge_write_realspace: unknown format "'//TRIM(fmt)//'"', &
                          calledby='wgauge_write_realspace')
      END SELECT
      CLOSE(iu)
   END SUBROUTINE wgauge_write_realspace

END MODULE m_wgauge_io
