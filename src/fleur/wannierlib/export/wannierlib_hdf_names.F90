!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> What the HDF files wannierlib writes are called, and which ones there are. One module
!> with nothing under it, so that everything which has to name one of them can reach it.
!>
!> WHY IT IS ITS OWN MODULE. The name belongs with the writer, the reader, and the pass
!> that stamps each file with the checkpoint it goes with -- and that last one cannot
!> reach into wannierlib_hdf_util, because hdf_util already reaches into it to write the
!> provenance attributes. Putting the name where hdf_util lives would make that a cycle,
!> so the spelling moved down here instead, below both. It depends on nothing but juDFT.
!>
!> THE ROSTER IS THE POINT, not just the spelling. wl_hdf_create refuses a name that is
!> not on it, so a fourth file cannot be added without appearing here -- and the stamp
!> pass walks this same list. Without that, a new exporter would write a file that comes
!> out with no provenance and no stamp, and nothing anywhere would say so: it would simply
!> not be in the list the stamp iterates. The guard turns a silent omission into a build
!> that runs once and aborts with the name in the message.
MODULE m_wannierlib_hdf_names
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wl_hdf_name, wl_hdf_check_what, WL_HDF_WHAT

   !> Every HDF file this writes, by the stem that goes after WF<n>_. Add one here first.
   CHARACTER(LEN=*), PARAMETER :: WL_HDF_WHAT(3) = &
      [CHARACTER(LEN=12) :: 'bmesh', 'operators_ks', 'modern_k']

CONTAINS

   !> WF<jspin>_<what>.hdf. The one place this is spelled, for the writer, the reader and
   !> the stamp alike: two spellings of one convention drift apart without anything
   !> noticing, and there were three of them before this module existed.
   FUNCTION wl_hdf_name(what, jspin) RESULT(filename)
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER,          INTENT(IN) :: jspin
      CHARACTER(LEN=80) :: filename
      WRITE (filename, '(a,i0,a)') 'WF', jspin, '_'//TRIM(what)//'.hdf'
   END FUNCTION wl_hdf_name

   !> Refuse a file that is not on the roster. Called where a file is CREATED, because
   !> that is the moment a new one would slip past the stamp pass unnoticed.
   SUBROUTINE wl_hdf_check_what(what, calledby)
      CHARACTER(LEN=*), INTENT(IN) :: what, calledby
      INTEGER :: i
      DO i = 1, SIZE(WL_HDF_WHAT)
         IF (TRIM(what) == TRIM(WL_HDF_WHAT(i))) RETURN
      END DO
      CALL juDFT_error("wannierlib: '"//TRIM(what)//"' is not one of the HDF files this "// &
                       "writes", hint='add it to WL_HDF_WHAT in wannierlib_hdf_names.F90, '// &
                       'which is also the list the checkpoint stamp walks', calledby=calledby)
   END SUBROUTINE wl_hdf_check_what

END MODULE m_wannierlib_hdf_names
