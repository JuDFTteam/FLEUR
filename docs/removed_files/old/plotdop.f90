!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_plotdop

use m_juDFT
use m_types

!+++++++++++++++++++++++++++++++++++++++++++++++++
!  plot the charge density for fleur  code output
!
!  if twodim = .false. a 3-D plot with nplot layers in z-direction
!  is constructed; the 3x3 matrix gives the 3 vectors of the cell ..
!  .gustav
!
!  Changed the input/output for easier use.
!  This subroutine uses the file plot_inp for input.
!  The old plotin-file is still supported by the old subroutine at the
!  end of the module
!                     Juelich, 21.1.06 DW
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++

CONTAINS

SUBROUTINE plotdop(oneD,stars,vacuum,sphhar,atoms,&
                   input,sym,cell,sliceplot,noco,cdnfname)

   USE m_outcdn
   USE m_loddop
   USE m_xsf_io
   USE m_cdn_io
   USE m_constants

   IMPLICIT NONE

   TYPE(t_oneD),                INTENT(IN)    :: oneD

   TYPE(t_stars),               INTENT(IN)    :: stars
   TYPE(t_vacuum),              INTENT(IN)    :: vacuum
   TYPE(t_sphhar),              INTENT(IN)    :: sphhar
   TYPE(t_atoms),               INTENT(IN)    :: atoms
   TYPE(t_input),               INTENT(IN)    :: input
   TYPE(t_sym),                 INTENT(IN)    :: sym
   TYPE(t_cell),                INTENT(IN)    :: cell
   TYPE(t_sliceplot),           INTENT(IN)    :: sliceplot
   TYPE(t_noco),                INTENT(IN)    :: noco
   CHARACTER(len=10), OPTIONAL, INTENT(IN)    :: cdnfname

!  .. Local Scalars ..
   REAL          :: tec,qint,fermiEnergyTemp,phi0,angss
   INTEGER       :: i,j,ix,iy,iz,jsp,na,nplo,iv,iflag,nfile
   INTEGER       :: nplot,nt,jm,jspin,numInFiles,numOutFiles
   LOGICAL       :: twodim,oldform,newform,l_qfix
   LOGICAL       :: cartesian,xsf,unwind,polar

!  .. Local Arrays ..
   TYPE(t_potden), ALLOCATABLE :: den(:)
   REAL, ALLOCATABLE    :: xdnout(:)
   REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3),help(3),qssc(3)
   INTEGER :: grid(3)
   REAL    :: rhocc(atoms%jmtd)
   REAL    :: point(3)
   CHARACTER (len=10), ALLOCATABLE :: cdnFilenames(:)
   CHARACTER (len=15), ALLOCATABLE :: outFilenames(:)
   CHARACTER (len=30)              :: filename
   CHARACTER (len=7)               :: textline
end subroutine
END MODULE m_plotdop
