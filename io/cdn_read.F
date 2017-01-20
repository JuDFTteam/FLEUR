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
     <                     ello,evac,epar,bkpt,wk,n_bands,n_size)
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
      REAL,    INTENT (INOUT) :: bkpt(:) !bkpt(3)
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
     <                    nmat,nv,ello,evdu,epar,kveclo,
     <                    k1,k2,k3,bkpt,wk,nbands,eig,z)

      USE m_eig66_io, ONLY : read_eig

      IMPLICIT NONE
!
! Arguments ...
!
      INTEGER, INTENT (IN) :: eig_id,irank,isize,ikpt,nbasfcn
      INTEGER, INTENT (IN) :: nvd,jspd,jspin
      INTEGER, INTENT (IN) :: noccbd,n_start,n_end
      LOGICAL, INTENT (IN) :: l_ss,l_noco
      INTEGER, INTENT (OUT) :: nbands,nmat
      REAL,    INTENT (OUT) :: wk

      INTEGER, INTENT (INOUT) :: k1(:,:),k2(:,:),k3(:,:) !k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (OUT) :: nv(:),kveclo(:) !nv(jspd),kveclo(nlotot)
      REAL,    INTENT (OUT) :: bkpt(:),eig(:) !bkpt(3),eig(neigd)
      REAL,    INTENT (INOUT) :: ello(:,:,:),evdu(:,:) !ello(nlod,ntypd,jspd),evdu(2,jspd)
      REAL,    INTENT (INOUT) :: epar(0:,:,:) !epar(0:lmaxd,ntypd,jspd)

      CLASS(*), INTENT (OUT) :: z(:,:) !z(nbasfcn,noccbd) !can be real/complex
!
! Local variables ...
!
      INTEGER :: iv,j,isp
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
     <                bk=bkpt,wk=wk,neig=nbands)
        DO isp = jspd,1,-1
          CALL read_eig(eig_id,ikpt,isp,
     <                    nmat=nmat,nv=nv(isp),k1=k1(:,isp),
     <                    k2=k2(:,isp),k3=k3(:,isp),kveclo=kveclo)
!         write(*,*) kveclo
        ENDDO
        CALL read_eig(
     >                 eig_id,ikpt,1,n_start=n_start,n_end=n_end,
     <                 eig=eig,nmat=nmat,z=z)
!
! For Non-Collinear, but no Spin-Spirals
!
      ELSEIF (l_noco) THEN
        CALL read_eig(
     >                eig_id,ikpt,1,
     <                bk=bkpt,wk=wk,neig=nbands,
     <                  nmat=nmat,nv=nv(1),k1=k1(:,1),
     <                  k2=k2(:,1),k3=k3(:,1),kveclo=kveclo)
        CALL read_eig(
     >                 eig_id,ikpt,1,n_start=n_start,n_end=n_end,
     <                 eig=eig,nmat=nmat,z=z)
         nv(jspd) = nv(1)
         DO iv = 1,nv(1)
            k1(iv,jspd) = k1(iv,1)
            k2(iv,jspd) = k2(iv,1)
            k3(iv,jspd) = k3(iv,1)
         ENDDO
!
! For Collinear
!
      ELSE
         CALL read_eig(eig_id,ikpt,jspin,
     <                   nmat=nmat)

         SELECT TYPE(z)
           TYPE IS (real)
             z=0
           TYPE IS (complex)
              z=0
          END SELECT
         CALL read_eig(
     >                  eig_id,ikpt,jspin,n_start=n_start,n_end=n_end,
     <                  bk=bkpt,wk=wk,
     <                 neig=nbands,eig=eig,
     <                   nv=nv(jspin),k1=k1(:,jspin),
     <                   k2=k2(:,jspin),k3=k3(:,jspin),kveclo=kveclo,
     <                  z=z(:,:noccbd))

      ENDIF
#ifdef CPP_MPI
!     The barrier here is not algorithmically required.
!     Nevertheless it improves the performance, presumably because it
!     synchronizes the processes and prevent the latency accumulating.
!     U.Alekseeva 2015

!Note: this does not work, two problems:
! a) MPI_COMM_WORLD should not be used 
! b) Not all processes might reach this point due to k-point parallelization
!      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr) 
#endif
      CALL timestop("cdn_read")
      
!      IF (nbands>neigd)   CALL juDFT_error("nbands.GT.neigd",calledby
!     +     ="cdn_read")
      IF (nv(jspin)>nvd)  CALL juDFT_error("nv.GT.nvd",calledby
     +     ="cdn_read")
      IF (nmat>nbasfcn)   CALL juDFT_error("nmat.GT.nbasfcn",calledby
     +     ="cdn_read")

      END SUBROUTINE cdn_read
      END MODULE m_cdnread
