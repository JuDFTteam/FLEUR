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
     <                     ello,evac,epar,wk,n_bands,n_size)
      USE m_eig66_io, ONLY : read_eig
      IMPLICIT NONE
!
! Arguments ...
!
      INTEGER, INTENT (IN) :: eig_id,irank,isize
      INTEGER, INTENT (IN) :: jspin,jspins
      LOGICAL, INTENT (IN) :: l_noco


      INTEGER, INTENT (OUT) :: n_size
      REAL,    INTENT (INOUT) :: wk

      INTEGER, INTENT (OUT) :: n_bands(0:) !n_bands(0:neigd)
      REAL,    INTENT (INOUT) :: ello(:,:,:),evac(:,:) !ello(nlod,ntypd,jspd),evac(2,jspd)
      REAL,    INTENT (INOUT) :: epar(0:,:,:) !epar(0:lmaxd,ntypd,jspd)

      INTEGER isp

      IF (l_noco) THEN
         CALL read_eig(eig_id,1,1,
     <                 neig=n_bands(1))
         DO isp = 1, jspins
           CALL read_eig(eig_id,1,isp,
     <             el=epar(:,:,isp),ello=ello(:,:,isp),evac=evac(:,isp))
         ENDDO

      ELSE
         CALL read_eig(eig_id,1,jspin,
     <                 neig=n_bands(1),el=epar(:,:,jspin),
     <                 ello=ello(:,:,jspin),evac=evac(:,jspin))

      ENDIF
c
c n_size is the number of records per k-point,
c n_bands(i) the number of ev's processed on n_rank=0...i-1
c
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
     <                    ello,evdu,epar,
     <                    wk,nbands,eig,zmat)

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
      REAL,    INTENT (OUT) :: wk

      REAL,    INTENT (OUT) :: eig(:) !bkpt(3),eig(neigd)
      REAL,    INTENT (INOUT) :: ello(:,:,:),evdu(:,:) !ello(nlod,ntypd,jspd),evdu(2,jspd)
      REAL,    INTENT (INOUT) :: epar(0:,:,:) !epar(0:lmaxd,ntypd,jspd)

      TYPE(t_zmat), INTENT (INOUT) :: zmat !z(nbasfcn,noccbd) !can be real/complex
!
! Local variables ...
!
      INTEGER :: iv,j,isp,nmat
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER mpiierr
#endif

!
! For Spin-Spirals
!
      CALL timestart("cdn_read")
      IF (l_ss) THEN
        CALL read_eig(eig_id,ikpt,1,
     <                wk=wk,neig=nbands)
        DO isp = jspd,1,-1
          CALL read_eig(eig_id,ikpt,isp,
     <                    nmat=nmat)
!         write(*,*) kveclo
        ENDDO
        CALL read_eig(
     >                 eig_id,ikpt,1,n_start=n_start,n_end=n_end,
     <                 eig=eig,nmat=nmat,zmat=zmat)
!
! For Non-Collinear, but no Spin-Spirals
!
      ELSEIF (l_noco) THEN
        CALL read_eig(
     >                eig_id,ikpt,1,
     <                wk=wk,neig=nbands,
     <                  nmat=nmat)
        CALL read_eig(
     >                 eig_id,ikpt,1,n_start=n_start,n_end=n_end,
     <                 eig=eig,nmat=nmat,zmat=zmat)
  !
! For Collinear
!
      ELSE
         CALL read_eig(eig_id,ikpt,jspin,
     <                   nmat=nmat)

         IF (zmat%l_real) THEN
            zmat%z_r=0
         ELSE
            zmat%z_c=0
         ENDIF
 
     
      CALL read_eig(
     >                  eig_id,ikpt,jspin,n_start=n_start,n_end=n_end,
     <                  wk=wk,
     <                  neig=nbands,eig=eig,
     <                  zmat=zmat)

      ENDIF
      CALL timestop("cdn_read")
      
!      IF (nbands>neigd)   CALL juDFT_error("nbands.GT.neigd",calledby
!     +     ="cdn_read")
      IF (nmat>nbasfcn)   CALL juDFT_error("nmat.GT.nbasfcn",calledby
     +     ="cdn_read")

      END SUBROUTINE cdn_read
      END MODULE m_cdnread
