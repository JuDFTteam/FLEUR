!TODO: Constants with capital letter
!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Contains required constants.
!
!> @author
!> Christian-Roman Gerhorst
!>
!> @brief
!> This module contains constant scalars and matrices which are needed within the program.
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpConstants.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_JPConstants

  use m_constants

  implicit none

  real,    parameter                   :: pi  =  4. * atan(1d0) !< The mathematical constant pi
  real,    parameter                   :: tpi =  8. * atan(1d0) !< Two times pi.
  real,    parameter                   :: fpi = 16. * atan(1d0) !< Four times pi.
  complex, parameter                   :: iu = cmplx(0, 1)      !< The imaginary unit.

  !Matrix syntax idea from http://stackoverflow.com/questions/3708307/how-to-initialize-two-dimensional-arrays-in-fortran
  complex, parameter, dimension(3, 3)  :: Tmatrix = transpose(reshape([ &
                                                                   & cmplx(1 / sqrt(2.), 0), cmplx(0, 0), cmplx(-1 / sqrt(2.), 0), &
                                                                   & -iu / sqrt(2.), cmplx(0, 0), -iu / sqrt(2.), &
                                                                   & cmplx(0, 0), cmplx(1, 0), cmplx(0, 0) &
                                                                   & ], [3, 3] )) !< Klüppelberg PhD thesis 4.18

  !Matrix syntax idea from http://stackoverflow.com/questions/3708307/how-to-initialize-two-dimensional-arrays-in-fortran
  complex, parameter, dimension(3, 3)  :: Tmatrix_transposed = reshape([ &
                                                                   & cmplx(1 / sqrt(2.), 0), cmplx(0, 0), cmplx(-1 / sqrt(2.), 0), &
                                                                   & -iu / sqrt(2.), cmplx(0, 0), -iu / sqrt(2.), &
                                                                   & cmplx(0, 0), cmplx(1, 0), cmplx(0, 0) &
                                                                   & ], [3, 3] ) !< Klüppelberg PhD thesis 4.18

  !Matrix syntax idea from http://stackoverflow.com/questions/3708307/how-to-initialize-two-dimensional-arrays-in-fortran
  complex, parameter, dimension(3, 3)  :: c_im = transpose(reshape([ &
                                                           & cmplx(sqrt(2 * pi / 3), 0), cmplx(0, 0), cmplx(-sqrt(2 * pi / 3), 0), &
                                                           & cmplx(0, sqrt(2 * pi / 3)), cmplx(0, 0), cmplx(0, sqrt(2 * pi / 3)), &
                                                           & cmplx(0, 0), cmplx(sqrt(4 * pi / 3), 0), cmplx(0, 0) &
                                                           & ], [3, 3] )) !< Klüppelberg PhD thesis 4.28

  !Matrix syntax idea from http://stackoverflow.com/questions/3708307/how-to-initialize-two-dimensional-arrays-in-fortran
  complex, parameter, dimension(3, 3)  :: c_mi = reshape([ &
                                                           & cmplx(sqrt(2 * pi / 3), 0), cmplx(0, 0), cmplx(-sqrt(2 * pi / 3), 0), &
                                                           & cmplx(0, sqrt(2 * pi / 3)), cmplx(0, 0), cmplx(0, sqrt(2 * pi / 3)), &
                                                           & cmplx(0, 0), cmplx(sqrt(4 * pi / 3), 0), cmplx(0, 0) &
                                                           & ], [3, 3] ) !< Klüppelberg PhD thesis 4.28

  ! Idea taken from Introduction to Modern Fortran. A lecture held by NickMacLaren in November 2007
  integer, parameter, dimension(3)    :: e1 = [1, 0, 0]
  integer, parameter, dimension(3)    :: e2 = [0, 1, 0]
  integer, parameter, dimension(3)    :: e3 = [0, 0, 1]
  integer, parameter, dimension(3, 3) :: idM = reshape([e1, e2, e3], [3, 3])

  ! A switch that activates/deactivates all Elk-Comparison output and related changes.
  logical, parameter                  :: compPhon = .false.!.true.

  ! A switch that enables the changes/fixes made by AN
  logical, parameter                  :: anfix = .false.!.true.

end module m_JPConstants
