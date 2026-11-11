!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Build the 2N Hamiltonian and position operator a collinear run owes orbitrans.
!>
!>  A collinear run wannierises each spin channel on its own, so the operator layer
!>  renders two 9x9 objects, WF1_* and WF2_*, while orbitrans reads one 18x18. That
!>  combination is physics, not a rename, so it cannot live in the renderer: the same
!>  rendering also feeds JuWIP and the standalone Wannier90 files, and those want the
!>  channels apart. It is assembled here instead, in the orbitrans consumer itself,
!>  from the files the renderer has already written.
!>
!>    hoppings(R) = diag( H_up(R), H_down(R) ) + Htr * SOC(R)
!>    position(R) = diag( A_up(R), A_down(R) )
!>
!>  in the order [up x num_wann, down x num_wann], which is the one the diagonal of
!>  S_z(R=0) shows and the one the spin indices of rssocmat.1 follow.
!>
!>  The spin-orbit term belongs INSIDE hoppings because orbitrans subtracts rather
!>  than adds: read_spin_orbit ends with H_nosoc = H_tot - H_soc, so it takes the
!>  hopping file for the full Hamiltonian and uses spin_orbit only to split the two
!>  apart. Handing it a spin-diagonal H instead leaves no matrix element between
!>  states of equal spin and sigma^Sx and sigma^Sy come out identically zero, with no
!>  error and with the orbital channel still looking healthy.
!>
!>  rssocmat.1 is in Hartree and the hopping file in eV, hence the conversion. Its
!>  orbital indices already run over all 2N functions and the (s1,s2) pair sits on
!>  top of them, so the element that counts is the one whose spin labels match those
!>  of the two functions. Treating the pair as a block offset instead leaves the
!>  up-down block untouched, which is the same silent zero again.
MODULE m_wannierlib_export_orbitrans_2n
   USE m_juDFT
   USE m_constants, ONLY: oUnit, hartree_to_ev_const
   USE m_wannierlib_export_paths, ONLY: wannierlib_export_file
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: orbitrans_assemble_2n

CONTAINS

   !> Replace the channel-one files in the export directory by their 2N assembly.
   !> Rank 0 only, and only once the last channel has been rendered.
   SUBROUTINE orbitrans_assemble_2n(num_wann)
      INTEGER, INTENT(IN) :: num_wann          !< functions per spin channel

      COMPLEX, ALLOCATABLE :: oup(:, :, :, :), odn(:, :, :, :), o2n(:, :, :, :)
      INTEGER, ALLOCATABLE :: rv(:, :), rv2(:, :)
      INTEGER :: nrpts, nrpts2, nw2
      LOGICAL :: l_soc, l_h, l_a

      nw2 = 2*num_wann

      !> Only what the run actually rendered. The consumer flag does not decide the
      !> operator list -- <operatorsRs> does -- so a run that asked for the Hamiltonian
      !> and not for the position leaves nothing to stack for the second one. Absent
      !> means not asked for; present on one channel and missing on the other is a
      !> broken run and says so below.
      INQUIRE (file=wannierlib_export_file('WF1_hr.dat'), exist=l_h)
      INQUIRE (file=wannierlib_export_file('WF1_r.dat'), exist=l_a)
      IF (.NOT. (l_h .OR. l_a)) RETURN

      !> --- H(R). The channel-one file already carries the orbitrans name.
      IF (l_h) THEN
         CALL orbitrans_read_op(wannierlib_export_file('WF1_hr.dat'), num_wann, 1, .FALSE., .FALSE., &
                                nrpts, rv, oup)
         CALL orbitrans_read_op(wannierlib_export_file('WF2_hr.dat'), num_wann, 1, .TRUE., .TRUE., &
                                nrpts2, rv2, odn)
         CALL orbitrans_same_mesh(nrpts, nrpts2, rv, rv2, 'WF2_hr.dat')
         CALL orbitrans_stack(oup, odn, num_wann, nrpts, 1, o2n)

         INQUIRE (file=wannierlib_export_file('rssocmat.1'), exist=l_soc)
         IF (l_soc) THEN
            CALL orbitrans_add_soc(o2n, num_wann, nrpts, rv)
         ELSE
            WRITE (oUnit, '(a)') 'wannierlib: WARNING -- no spin-orbit operator to fold into '// &
               'hoppings. The Hamiltonian stays diagonal in spin and the transverse spin '// &
               'conductivities will come out identically zero.'
         END IF
         CALL orbitrans_write_op(wannierlib_export_file('WF1_hr.dat'), nw2, nrpts, 1, rv, o2n)
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: assembled hoppings as ', nw2, 'x', nw2, &
            ' from the two spin channels'
      END IF

      !> --- A(R), the same stacking and no extra term. Its R mesh is checked against the
      !> file itself rather than against the Hamiltonian's, which a run without one has not
      !> read.
      IF (l_a) THEN
         CALL orbitrans_read_op(wannierlib_export_file('WF1_r.dat'), num_wann, 3, .FALSE., .FALSE., &
                                nrpts, rv, oup)
         CALL orbitrans_read_op(wannierlib_export_file('WF2_r.dat'), num_wann, 3, .TRUE., .FALSE., &
                                nrpts2, rv2, odn)
         CALL orbitrans_same_mesh(nrpts, nrpts2, rv, rv2, 'WF2_r.dat')
         CALL orbitrans_stack(oup, odn, num_wann, nrpts, 3, o2n)
         CALL orbitrans_write_op(wannierlib_export_file('WF1_r.dat'), nw2, nrpts, 3, rv, o2n)
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: assembled position as ', nw2, 'x', nw2, &
            ' from the two spin channels'
      END IF
   END SUBROUTINE orbitrans_assemble_2n

   !> One real-space operator file into (nw,nw,nrpts,ncomp), with its R vectors.
   !>
   !> The channel-one files have had their header taken off and their degeneracies moved
   !> to `ndegen`, because that is what the consumer reads; WF2_* keeps the Wannier90
   !> header, so the number of R vectors is read from it rather than guessed, and the
   !> degeneracy block that follows is ceil(nrpts/15) lines wide.
   SUBROUTINE orbitrans_read_op(fname, nw, ncomp, l_header, l_ndegen, nrpts, rv, o)
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: nw, ncomp
      LOGICAL, INTENT(IN) :: l_header    !< the Wannier90 three-line header is still on
      LOGICAL, INTENT(IN) :: l_ndegen    !< and behind it the degeneracy block (hr only)
      INTEGER, INTENT(OUT) :: nrpts
      INTEGER, ALLOCATABLE, INTENT(OUT) :: rv(:, :)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: o(:, :, :, :)

      INTEGER :: iu, ios, irpt, i, j, c, k, nw_file, ndlines, ir(3), ii, jj
      REAL :: val(2*ncomp)
      CHARACTER(LEN=256) :: line

      OPEN (newunit=iu, file=TRIM(fname), status='old', action='read', iostat=ios)
      IF (ios /= 0) CALL juDFT_error('<export orbitrans="T"> cannot open '//TRIM(fname)// &
                                     ': a collinear run has to leave both spin channels behind', &
                                     calledby='orbitrans_read_op')

      IF (l_header) THEN
         READ (iu, '(a)') line
         READ (iu, *) nw_file
         READ (iu, *) nrpts
         IF (nw_file /= nw) CALL juDFT_error('channel two of '//TRIM(fname)//' has a '// &
                                             'different number of Wannier functions', calledby='orbitrans_read_op')
         IF (l_ndegen) THEN
            ndlines = (nrpts + 14)/15
            DO k = 1, ndlines
               READ (iu, '(a)') line
            END DO
         END IF
      ELSE
         CALL orbitrans_count_rpts(TRIM(fname), nw, nrpts)
      END IF

      ALLOCATE (o(nw, nw, nrpts, ncomp), rv(3, nrpts))
      DO irpt = 1, nrpts
         DO j = 1, nw
            DO i = 1, nw
               READ (iu, *, iostat=ios) ir, ii, jj, val
               IF (ios /= 0) CALL juDFT_error(TRIM(fname)//' ended early: it does not hold '// &
                                              'nrpts*nw*nw records', calledby='orbitrans_read_op')
               IF (ii /= i .OR. jj /= j) &
                  CALL juDFT_error(TRIM(fname)//' is not in the order this reader assumes '// &
                                   '(R outermost, then column, then row)', calledby='orbitrans_read_op')
               IF (i == 1 .AND. j == 1) rv(:, irpt) = ir
               DO c = 1, ncomp
                  o(i, j, irpt, c) = CMPLX(val(2*c - 1), val(2*c))
               END DO
            END DO
         END DO
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_read_op

   !> How many R vectors a headerless file holds, from its record count.
   SUBROUTINE orbitrans_count_rpts(fname, nw, nrpts)
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: nw
      INTEGER, INTENT(OUT) :: nrpts

      INTEGER :: iu, ios, n
      CHARACTER(LEN=256) :: line

      n = 0
      OPEN (newunit=iu, file=TRIM(fname), status='old', action='read')
      DO
         READ (iu, '(a)', iostat=ios) line
         IF (ios /= 0) EXIT
         IF (LEN_TRIM(line) > 0) n = n + 1
      END DO
      CLOSE (iu)
      IF (MOD(n, nw*nw) /= 0) CALL juDFT_error(TRIM(fname)//' does not hold a whole number '// &
                                               'of nw x nw blocks', calledby='orbitrans_count_rpts')
      nrpts = n/(nw*nw)
   END SUBROUTINE orbitrans_count_rpts

   !> Both channels have to live on the same Wigner-Seitz mesh, or the blocks being put
   !> on one diagonal belong to different R vectors. Same mesh, same order, every vector.
   SUBROUTINE orbitrans_same_mesh(nrpts, nrpts2, rv, rv2, who)
      INTEGER, INTENT(IN) :: nrpts, nrpts2, rv(:, :), rv2(:, :)
      CHARACTER(LEN=*), INTENT(IN) :: who

      INTEGER :: irpt

      IF (nrpts /= nrpts2) CALL juDFT_error('the two spin channels disagree on the number '// &
                                            'of R vectors ('//TRIM(who)//')', calledby='orbitrans_same_mesh')
      DO irpt = 1, nrpts
         IF (ANY(rv(:, irpt) /= rv2(:, irpt))) &
            CALL juDFT_error('the two spin channels list the R vectors in a different order ('// &
                             TRIM(who)//')', calledby='orbitrans_same_mesh')
      END DO
   END SUBROUTINE orbitrans_same_mesh

   !> diag(up, down) in the order [up x nw, down x nw]. The off-diagonal blocks are the
   !> ones the spin-orbit term later fills; without it they stay zero, which is what a
   !> Hamiltonian that does not mix the channels means.
   SUBROUTINE orbitrans_stack(oup, odn, nw, nrpts, ncomp, o2n)
      COMPLEX, INTENT(IN) :: oup(:, :, :, :), odn(:, :, :, :)
      INTEGER, INTENT(IN) :: nw, nrpts, ncomp
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: o2n(:, :, :, :)

      IF (ALLOCATED(o2n)) DEALLOCATE (o2n)
      ALLOCATE (o2n(2*nw, 2*nw, nrpts, ncomp), source=CMPLX(0.0, 0.0))
      o2n(1:nw, 1:nw, :, :) = oup
      o2n(nw + 1:2*nw, nw + 1:2*nw, :, :) = odn
   END SUBROUTINE orbitrans_stack

   !> Fold rssocmat.1 into the stacked Hamiltonian, in eV.
   !>
   !> The file is nrpts x 2N x 2N x 2 x 2: the orbital indices already run over the whole
   !> 2N basis and carry the spinor pair on top. Function i belongs to the up channel when
   !> i <= nw, so the element that is physically its own is the one with s1 and s2 equal to
   !> those labels; the other three of each quartet belong to a different pair of functions.
   SUBROUTINE orbitrans_add_soc(h2n, nw, nrpts, rv)
      COMPLEX, INTENT(INOUT) :: h2n(:, :, :, :)
      INTEGER, INTENT(IN) :: nw, nrpts, rv(:, :)

      INTEGER :: iu, ios, irpt, k, nw2, ir(3), i, j, s1, s2, si, sj
      REAL :: re, im, cross
      CHARACTER(LEN=:), ALLOCATABLE :: fname

      nw2 = 2*nw
      fname = wannierlib_export_file('rssocmat.1')
      OPEN (newunit=iu, file=fname, status='old', action='read', iostat=ios)
      IF (ios /= 0) CALL juDFT_error('cannot open '//fname, calledby='orbitrans_add_soc')

      DO irpt = 1, nrpts
         DO k = 1, nw2*nw2*4
            READ (iu, *, iostat=ios) ir, i, j, s1, s2, re, im
            IF (ios /= 0) CALL juDFT_error(fname//' ended early: it does not hold '// &
                                           'nrpts*2N*2N*4 records', calledby='orbitrans_add_soc')
            IF (ANY(ir /= rv(:, irpt))) &
               CALL juDFT_error(fname//' runs over the R vectors in a different order than '// &
                                'the hopping file', calledby='orbitrans_add_soc')
            si = MERGE(1, 2, i <= nw)
            sj = MERGE(1, 2, j <= nw)
            IF (s1 /= si .OR. s2 /= sj) CYCLE
            h2n(i, j, irpt, 1) = h2n(i, j, irpt, 1) + hartree_to_ev_const*CMPLX(re, im)
         END DO
      END DO
      CLOSE (iu)

      !> The whole point of the fold is the block that connects the two channels. If it is
      !> still empty the spin indices did not line up, and the run would end in a transverse
      !> spin conductivity of exactly zero with nothing to show for it.
      cross = 0.0
      DO irpt = 1, nrpts
         DO j = 1, nw2
            DO i = 1, nw2
               IF ((i <= nw) .NEQV. (j <= nw)) cross = MAX(cross, ABS(h2n(i, j, irpt, 1)))
            END DO
         END DO
      END DO
      IF (cross <= 0.0) CALL juDFT_error('the spin-orbit term left the up-down block of '// &
                                         'hoppings empty, so it did not enter the Hamiltonian', &
                                         calledby='orbitrans_add_soc')
      WRITE (oUnit, '(a,f12.6,a)') 'wannierlib: spin-orbit folded into hoppings, '// &
         'largest up-down element ', cross, ' eV'
   END SUBROUTINE orbitrans_add_soc

   !> Write it back in the layout the file came in: R outermost, then column, then row,
   !> and no header, which is what the consumer reads from the first line.
   SUBROUTINE orbitrans_write_op(fname, nw, nrpts, ncomp, rv, o)
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: nw, nrpts, ncomp, rv(:, :)
      COMPLEX, INTENT(IN) :: o(:, :, :, :)

      INTEGER :: iu, irpt, i, j, c
      REAL :: val(2*ncomp)

      OPEN (newunit=iu, file=TRIM(fname), status='replace')
      DO irpt = 1, nrpts
         DO j = 1, nw
            DO i = 1, nw
               DO c = 1, ncomp
                  val(2*c - 1) = REAL(o(i, j, irpt, c))
                  val(2*c) = AIMAG(o(i, j, irpt, c))
               END DO
               WRITE (iu, '(5i5,6f12.6)') rv(:, irpt), i, j, val
            END DO
         END DO
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_op

END MODULE m_wannierlib_export_orbitrans_2n
