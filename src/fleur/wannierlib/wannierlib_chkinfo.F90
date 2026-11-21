!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> What the checkpoint cannot say about itself, carried by the files that are read
!> beside it.
!>
!> WHY IT IS NEEDED. Wannier90 checks ten fields when it loads a .chk -- num_bands,
!> num_wann, num_kpts, nntot, mp_grid, kpt_latt, the two lattices, exclude_bands and
!> num_exclude_bands -- and that answers "is this checkpoint shaped like this case". It
!> does not answer "is it from this run". The energy windows never enter the file at all;
!> only their consequence does, the lwindow mask, which is read and imposed. So a .chk
!> from a run with different disWinMin, different disFrozMax, or with minBand and maxBand
!> shifted in step -- numBands unchanged -- passes all ten checks and hands over a gauge
!> built on other bands, with no error anywhere.
!>
!> WHO CHECKS. Not FLEUR. The code that has the problem is whatever reads the checkpoint
!> next -- orbitrans, a script, another code -- and it is the one holding both artefacts
!> at once. FLEUR's job is to state, beside the gauge it hands out, what that gauge was
!> built from, and to say which checkpoint it goes with. Refusing on FLEUR's side would
!> only cover FLEUR's own restart and would leave every outside consumer where it was.
!>
!> WHERE IT LIVES. In the root group of WF<n>_bmesh.hdf, WF<n>_operators_ks.hdf and
!> WF<n>_modern_k.hdf, not in a file of its own. A separate .chkinfo would have to be
!> paired with its .chk by file name, which is the same kind of trust the stamp exists to
!> replace; and the files that need the information are exactly these, because they are
!> what a consumer reads the checkpoint alongside. The .chk itself is not an option:
!> Wannier90 owns that format and we read it, we do not write to it.
!>
!> HOW IT READS. Ten integer and real attributes on the root group, plus two character
!> ones. FLEUR's attribute writer stores a string as an ARRAY OF ONE-CHARACTER STRINGS,
!> which is what every character attribute in every FLEUR HDF file looks like: h5dump
!> prints it as a word, but h5py hands back array([b'W', b'F', b'1', ...]) and comparing
!> that to a string is quietly False rather than an error. Joining is the whole of it --
!> b"".join(handle.attrs["route"]).decode() -- and it is said here because the place a
!> consumer looks for the convention is the file that writes it.
!>
!> THE STAMP is CRC32 over the bytes of the .chk, the reflected-polynomial variety zlib
!> computes, so checking it costs a consumer one line -- zlib.crc32(open(p,'rb').read())
!> in Python -- and no agreement with us about anything but the number. It is what ties
!> the two files together when every field they share already agrees, which is precisely
!> the shifted-band case above.
MODULE m_wannierlib_chkinfo
   USE m_juDFT
   USE m_wannierlib_hdf_names, ONLY: wl_hdf_name, WL_HDF_WHAT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wl_chkinfo_set, wl_chkinfo_stamp, wl_crc32_file
#ifdef CPP_HDF
   PUBLIC :: wl_chkinfo_write
#endif

   TYPE :: t_wl_chkinfo
      LOGICAL           :: l_set        = .FALSE.
      CHARACTER(LEN=16) :: route        = 'unset'
      INTEGER           :: num_bands    = 0
      INTEGER           :: min_band     = 0
      INTEGER           :: max_band     = 0
      INTEGER           :: num_wann     = 0
      INTEGER           :: nranks       = 0
      REAL              :: dis_win_min  = 0.0
      REAL              :: dis_win_max  = 0.0
      REAL              :: dis_froz_min = 0.0
      REAL              :: dis_froz_max = 0.0
   END TYPE t_wl_chkinfo

   !> Module state on purpose. Every one of these is settled once, before the first export
   !> runs, and is read again three layers down inside wl_hdf_root, which has no business
   !> growing ten arguments that none of its callers care about. Set once per channel.
   TYPE(t_wl_chkinfo), SAVE :: info

CONTAINS

   !> Called once per spin channel, after the windows have been resolved and before
   !> anything is written. The values are the EFFECTIVE ones: a window left out of the
   !> input is filled in from the bands themselves, and what belongs here is what the run
   !> used, not what the user typed.
   SUBROUTINE wl_chkinfo_set(route, num_bands, min_band, max_band, num_wann, nranks, &
                             dis_win_min, dis_win_max, dis_froz_min, dis_froz_max)
      CHARACTER(LEN=*), INTENT(IN) :: route
      INTEGER,          INTENT(IN) :: num_bands, min_band, max_band, num_wann, nranks
      REAL,             INTENT(IN) :: dis_win_min, dis_win_max, dis_froz_min, dis_froz_max
      info%l_set = .TRUE.
      info%route = route
      info%num_bands = num_bands
      info%min_band = min_band
      info%max_band = max_band
      info%num_wann = num_wann
      info%nranks = nranks
      info%dis_win_min = dis_win_min
      info%dis_win_max = dis_win_max
      info%dis_froz_min = dis_froz_min
      info%dis_froz_max = dis_froz_max
   END SUBROUTINE wl_chkinfo_set

   !> CRC32 of a file, the reflected-polynomial form (0xEDB88320) that zlib, gzip and
   !> PNG use. Zero if the file is not there, which is also what a run that handed out no
   !> checkpoint writes -- the two cases mean the same thing to a reader: nothing to pair
   !> with. ISHFT is a logical shift in Fortran, so the sign bit of a default INTEGER
   !> carries no meaning here and the table algorithm is the textbook one.
   INTEGER FUNCTION wl_crc32_file(path) RESULT(crc)
      CHARACTER(LEN=*), INTENT(IN) :: path
      INTEGER, PARAMETER :: POLY = -306674912          ! 0xEDB88320 as a signed 32-bit int
      INTEGER, PARAMETER :: CHUNK = 65536
      INTEGER            :: tbl(0:255)
      INTEGER            :: i, j, c, iu, ios, n, b
      INTEGER(KIND=1)    :: buf(CHUNK)
      INTEGER(KIND=8)    :: fsize, pos
      LOGICAL            :: l_ex

      crc = 0
      fsize = -1
      INQUIRE (file=path, exist=l_ex, size=fsize)
      IF (.NOT. l_ex .OR. fsize <= 0) RETURN

      DO i = 0, 255
         c = i
         DO j = 1, 8
            IF (IAND(c, 1) /= 0) THEN
               c = IEOR(ISHFT(c, -1), POLY)
            ELSE
               c = ISHFT(c, -1)
            END IF
         END DO
         tbl(i) = c
      END DO

      OPEN (NEWUNIT=iu, FILE=path, ACCESS='stream', FORM='unformatted', STATUS='old', &
            ACTION='read', IOSTAT=ios)
      IF (ios /= 0) RETURN

      crc = NOT(0)
      pos = 1
      DO WHILE (pos <= fsize)
         n = INT(MIN(INT(CHUNK, 8), fsize - pos + 1))
         READ (iu, POS=pos, IOSTAT=ios) buf(1:n)
         IF (ios /= 0) EXIT
         DO i = 1, n
            b = IAND(INT(buf(i)), 255)
            crc = IEOR(ISHFT(crc, -8), tbl(IAND(IEOR(crc, b), 255)))
         END DO
         pos = pos + n
      END DO
      CLOSE (iu)
      crc = NOT(crc)
   END FUNCTION wl_crc32_file

#ifdef CPP_HDF
   !> The parameters, onto a root group that is already open. Called from wl_hdf_root, so
   !> every file we write carries them and none of the exporters has to remember to.
   SUBROUTINE wl_chkinfo_write(gid)
      USE hdf5
      USE m_hdf_tools
      INTEGER(HID_T), INTENT(IN) :: gid
      IF (.NOT. info%l_set) RETURN
      CALL io_write_attchar0(gid, 'route', TRIM(info%route))
      CALL io_write_attint0(gid, 'numBands', info%num_bands)
      CALL io_write_attint0(gid, 'minBand', info%min_band)
      CALL io_write_attint0(gid, 'maxBand', info%max_band)
      CALL io_write_attint0(gid, 'numWann', info%num_wann)
      CALL io_write_attint0(gid, 'nranks', info%nranks)
      !> Hartree and absolute, as FLEUR holds them -- NOT relative to the Fermi level, and
      !> not the eV the operator datasets are printed in.
      CALL io_write_attreal0(gid, 'disWinMin', info%dis_win_min)
      CALL io_write_attreal0(gid, 'disWinMax', info%dis_win_max)
      CALL io_write_attreal0(gid, 'disFrozMin', info%dis_froz_min)
      CALL io_write_attreal0(gid, 'disFrozMax', info%dis_froz_max)
   END SUBROUTINE wl_chkinfo_write
#endif

   !> Add the checkpoint's fingerprint to every file of this channel that was written.
   !>
   !> Separate from wl_chkinfo_write, and late, because the files are not written at one
   !> time: WF<n>_bmesh.hdf and WF<n>_operators_ks.hdf go out BEFORE run_w90, when there
   !> is no checkpoint to hash yet, and WF<n>_modern_k.hdf after it. One pass once the
   !> channel is finished is the only point at which all three exist, so it is also the
   !> only place this can be one code path instead of two.
   !>
   !> Rank 0 alone: it is the rank that wrote the files, and reopening an HDF5 file for
   !> write from several ranks at once is not something any of them would survive.
   !>
   !> The roster and the spelling both come from wannierlib_hdf_names, the same ones
   !> wl_hdf_create works off, so a file this walks past is a file that was never created.
   !> A file that is MISSING is the normal case: each is written only if its own <export>
   !> flag or operator block asked for it.
   SUBROUTINE wl_chkinfo_stamp(jspin, chkfile, irank)
#ifdef CPP_HDF
      USE hdf5
      USE m_hdf_tools
#endif
      USE m_constants, ONLY: oUnit
      INTEGER,          INTENT(IN) :: jspin
      CHARACTER(LEN=*), INTENT(IN) :: chkfile
      INTEGER,          INTENT(IN) :: irank
#ifdef CPP_HDF
      INTEGER(HID_T)    :: fid, gid
      CHARACTER(LEN=80) :: fn
      INTEGER           :: i, e, crc
      LOGICAL           :: l_ex

      IF (irank /= 0) RETURN
      crc = wl_crc32_file(chkfile)
      DO i = 1, SIZE(WL_HDF_WHAT)
         fn = wl_hdf_name(WL_HDF_WHAT(i), jspin)
         INQUIRE (file=TRIM(fn), exist=l_ex)
         IF (.NOT. l_ex) CYCLE
         CALL h5fopen_f(TRIM(fn), H5F_ACC_RDWR_F, fid, e)
         IF (e /= 0) CYCLE
         CALL h5gopen_f(fid, '/', gid, e)
         CALL io_write_attint0(gid, 'chkCRC32', crc)
         CALL io_write_attchar0(gid, 'chkFile', TRIM(chkfile))
         CALL h5gclose_f(gid, e)
         CALL h5fclose_f(fid, e)
         WRITE (oUnit, '(a,i0)') 'wannierlib: '//TRIM(fn)//' stamped with the checkpoint ' &
            //TRIM(chkfile)//', CRC32 ', crc
      END DO
#endif
   END SUBROUTINE wl_chkinfo_stamp

END MODULE m_wannierlib_chkinfo
