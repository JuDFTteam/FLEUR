!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_cdnread
      use m_juDFT
!------------------------------------------------------------------------
!     Two subroutines to read in the information stored on the 'eig'-file
!     cdn_read0: reads energy-parameters and number of bands stored/pe
!     cdn_read : reads finally the eigenvalues and vectors
!                                                          gb`04
!------------------------------------------------------------------------
      CONTAINS
      SUBROUTINE cdn_read0(
     >                     eig_id,
     >                     irank,isize,jspin,jspins,
     >                     l_noco,
     <                     n_bands,n_size)
      USE m_eig66_io, ONLY : read_eig
      IMPLICIT NONE

! Arguments ...
      INTEGER, INTENT (IN) :: eig_id,irank,isize
      INTEGER, INTENT (IN) :: jspin,jspins
      LOGICAL, INTENT (IN) :: l_noco


      INTEGER, INTENT (OUT) :: n_size
 
      INTEGER, INTENT (OUT) :: n_bands(0:) !n_bands(0:neigd)
 
      INTEGER isp

      isp = MERGE(1,jspin,l_noco)

      CALL read_eig(eig_id,1,isp,neig=n_bands(1))

c n_size is the number of records per k-point,
c n_bands(i) the number of ev's processed on n_rank=0...i-1

      n_size = 1

      n_bands(0) = 0

      END SUBROUTINE cdn_read0
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
      SUBROUTINE cdn_read(
     >                    eig_id,nvd,jspd,irank,isize,
     >                    ikpt,jspin,nbasfcn,l_ss,l_noco,
     >                    noccbd,n_start,n_end,
     <                    nbands,eig,zmat)

      USE m_eig66_io, ONLY : read_eig
      USE m_types
 
      IMPLICIT NONE
!
! Arguments ...
!
      INTEGER, INTENT (IN) :: eig_id,irank,isize,ikpt,nbasfcn
      INTEGER, INTENT (IN) :: nvd,jspd,jspin
      INTEGER, INTENT (IN) :: noccbd,n_start,n_end
      LOGICAL, INTENT (IN) :: l_ss,l_noco
      INTEGER, INTENT (OUT) :: nbands
   
      REAL,    INTENT (OUT) :: eig(:) !bkpt(3),eig(neigd)
   
      TYPE(t_mat), INTENT (INOUT) :: zmat !z(nbasfcn,noccbd) !can be real/complex

! Local variables ...
      INTEGER :: isp,i
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER mpiierr
#endif

      isp = MERGE(1,jspin,l_noco)

      CALL timestart("cdn_read")

      CALL read_eig(eig_id,ikpt,isp,list=[(i,i=n_start,n_end)],
     <              neig=nbands,eig=eig,zmat=zmat)

      CALL timestop("cdn_read")
      
!      IF (nbands>neigd)   CALL juDFT_error("nbands.GT.neigd",calledby
!     +     ="cdn_read")
    
      END SUBROUTINE cdn_read
      END MODULE m_cdnread
