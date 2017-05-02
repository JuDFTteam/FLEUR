!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eig66_io
#include "juDFT_env.h"
  USE m_types
  USE m_eig66_data
  IMPLICIT NONE
  PRIVATE

  PUBLIC open_eig,close_eig
  PUBLIC read_eig, write_eig
  PUBLIC read_dos,write_dos
CONTAINS

  FUNCTION open_eig(mpi_comm,nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,l_noco,l_create,l_real,l_soc,l_readonly,n_size,mode_in,filename,layers,nstars,ncored,nsld,nat,l_dos,l_mcd,l_orb)RESULT(id)
    USE m_eig66_hdf,ONLY:open_eig_hdf=>open_eig
    USE m_eig66_DA ,ONLY:open_eig_DA=>open_eig
    USE m_eig66_mem,ONLY:open_eig_mem=>open_eig
    USE m_eig66_MPI,ONLY:open_eig_mpi=>open_eig
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,mpi_comm
    LOGICAL,INTENT(IN)          :: l_noco,l_readonly,l_create,l_real,l_soc
    INTEGER,INTENT(IN),OPTIONAL :: n_size,mode_in
    LOGICAL,INTENT(IN),OPTIONAL :: l_dos,l_mcd,l_orb
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
    INTEGER,INTENT(IN),OPTIONAL :: layers,nstars,ncored,nsld,nat
    INTEGER:: id,mode

    INTEGER:: neig_local,isize,err
    if (l_soc) THEN
       neig_local=2*neig
    else
       neig_local=neig
    endif
    mode=-1
    IF (PRESENT(mode_in)) mode=mode_in

    IF (mode<0) THEN
       !Use default mode
#ifdef CPP_MPI
       mode=MPI_mode
#else
       mode=MEM_mode
#endif
       !check if default was given on command-line
       IF (juDFT_was_argument("-mpi")) mode=MPI_mode
       IF (juDFT_was_argument("-mem")) mode=MEM_mode
       IF (juDFT_was_argument("-da")) mode=DA_mode
       IF (juDFT_was_argument("-hdf")) mode=HDF_mode
    ENDIF
    !Check if mode is available
#ifndef CPP_MPI
    IF (mode==MPI_mode) CALL juDFT_error("MPI-mode not available. Recompile with CPP_MPI",calledby="eig66_io")
#else
    CALL MPI_COMM_SIZE(mpi_comm,isize,err)
    IF (isize>1.AND.((mode==DA_mode.OR.mode==mem_mode))) &
         CALL juDFT_error("In a parallel calculation MEM/DA-mode are not available",calledby="eig66_io")
#endif
#ifndef CPP_HDF
    IF (mode==HDF_mode) CALL juDFT_error("HDF-mode not available. Recompile with CPP_HDF",calledby="eig66_io")
#endif

    id=eig66_data_newid(mode)

    !PRINT *,"open_eig:",id,mode

    CALL timestart("Open file/memory for IO of eig66")
    SELECT CASE (eig66_data_mode(id))
    CASE (DA_mode)
       CALL open_eig_DA(id,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,nlotot,l_create,l_real,l_soc,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)
    CASE (hdf_mode)
       CALL open_eig_HDF(id,mpi_comm,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,l_create,l_real,l_soc,nlotot,l_readonly,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)
    CASE (mem_mode)
       CALL open_eig_MEM(id,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,l_create,l_real,l_soc,nlotot,l_noco,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)
    CASE (mpi_mode)
       CALL open_eig_MPI(id,mpi_comm,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,l_create,l_real,l_soc,nlotot,l_noco,n_size,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)
    CASE DEFAULT
       CALL juDFT_error("Invalid IO-mode in eig66_io")
    END SELECT
    CALL timestop("Open file/memory for IO of eig66")
  END FUNCTION open_eig

  SUBROUTINE close_eig(id,filename)
    USE m_eig66_hdf,ONLY:close_eig_hdf=>close_eig
    USE m_eig66_DA ,ONLY:close_eig_DA=>close_eig
    USE m_eig66_mem,ONLY:close_eig_MEM=>close_eig
    USE m_eig66_MPI,ONLY:close_eig_MPI=>close_eig
    IMPLICIT NONE
    INTEGER,INTENT(IN)                   :: id
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
    INTEGER  :: mode
    mode=eig66_data_mode(id)
    !PRINT*,"close_eig:",id,mode
    SELECT CASE (mode)
    CASE (DA_mode)
       CALL close_eig_DA(id,filename)
    CASE (hdf_mode)
       CALL close_eig_HDF(id,filename)
    CASE (mem_mode)
       CALL close_eig_Mem(id,filename=filename)
    CASE (MPI_mode)
       CALL close_eig_MPI(id,filename=filename)
    CASE (-1)
       CALL juDFT_error("ID not assigned in close_eig",calledby="eig66_io")
    END SELECT

  END SUBROUTINE close_eig

  SUBROUTINE read_eig(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,zmat)
    USE m_eig66_hdf,ONLY:read_eig_hdf=>read_eig
    USE m_eig66_DA ,ONLY:read_eig_DA=>read_eig
    USE m_eig66_mem,ONLY:read_eig_mem=>read_eig
    USE m_eig66_MPI,ONLY:read_eig_MPI=>read_eig
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: id,nk,jspin
    INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
    INTEGER, INTENT(OUT),OPTIONAL  :: neig
    REAL,    INTENT(OUT),OPTIONAL  :: eig(:)
    INTEGER, INTENT(OUT),OPTIONAL  :: k1(:),k2(:),k3(:),kveclo(:)
    REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
    REAL,    INTENT(OUT),OPTIONAL  :: bk(:),wk
    INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
    TYPE(t_zMAT),INTENT(INOUT),OPTIONAL  :: zmat
    INTEGER::n
    CALL timestart("IO (read)")
    SELECT CASE (eig66_data_mode(id))
    CASE (DA_mode)
       CALL read_eig_DA(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,zmat)
    CASE (hdf_mode)
       CALL read_eig_hdf(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,zmat)
    CASE (mem_mode)
       CALL read_eig_mem(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,zmat)
    CASE (mpi_mode)
       CALL read_eig_mpi(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,zmat)
    CASE (-1)
       CALL juDFT_error("Could not read eig-file before opening")
    END SELECT
    CALL timestop("IO (read)")
  END SUBROUTINE read_eig

  SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,zmat)
    USE m_eig66_hdf,ONLY:write_eig_hdf=>write_eig
    USE m_eig66_DA ,ONLY:write_eig_DA=>write_eig
    USE m_eig66_mem,ONLY:write_eig_MEM=>write_eig
    USE m_eig66_MPI,ONLY:write_eig_MPI=>write_eig
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id,nk,jspin
    REAL,    INTENT(IN),OPTIONAL :: wk
    INTEGER, INTENT(IN),OPTIONAL :: neig,neig_total,nv,nmat,nlotot,n_start,n_end
    INTEGER, INTENT(IN),OPTIONAL :: k1(:),k2(:),k3(:),kveclo(:)
    REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:),evac(2),ello(:,:)
    TYPE(t_zMat),INTENT(IN),OPTIONAL :: zmat
    CALL timestart("IO (write)")
    SELECT CASE (eig66_data_mode(id))
    CASE (da_mode)
       CALL write_eig_DA(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,zmat)
    CASE (hdf_mode)
       CALL write_eig_HDF(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,zmat)
    CASE (mem_mode)
       CALL write_eig_Mem(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,zmat)
    CASE (MPI_mode)
       CALL write_eig_MPI(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,zmat)
    CASE (-1)
       CALL juDFT_error("Could not write eig-file before opening")
    END SELECT
    CALL timestop("IO (write)")
  END SUBROUTINE write_eig

  SUBROUTINE write_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    USE m_eig66_hdf,ONLY:write_dos_hdf=>write_dos
    USE m_eig66_DA ,ONLY:write_dos_DA=>write_dos
    USE m_eig66_mem,ONLY:write_dos_MEM=>write_dos
    USE m_eig66_MPI,ONLY:write_dos_MPI=>write_dos
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(IN)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(IN)           :: qstars(:,:,:,:)
    INTEGER,INTENT(IN)           :: ksym(:),jsym(:)
    REAL,INTENT(IN),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(IN),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
    CALL timestart("IO (dos-write)")
    SELECT CASE (eig66_data_mode(id))
    CASE (da_mode)
       CALL write_dos_DA(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (hdf_mode)
       CALL write_dos_HDF(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (mem_mode)
       CALL write_dos_Mem(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (MPI_mode)
       CALL write_dos_MPI(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (-1)
       CALL juDFT_error("Could not write eig-file before opening")
    END SELECT
    CALL timestop("IO (dos-write)")
  END SUBROUTINE write_dos


  SUBROUTINE read_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    USE m_eig66_hdf,ONLY:read_dos_hdf=>read_dos
    USE m_eig66_DA ,ONLY:read_dos_DA=>read_dos
    USE m_eig66_mem,ONLY:read_dos_MEM=>read_dos
    USE m_eig66_MPI,ONLY:read_dos_MPI=>read_dos
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(OUT)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(OUT)           :: qstars(:,:,:,:)
    INTEGER,INTENT(OUT)           :: ksym(:),jsym(:)
    REAL,INTENT(OUT),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(OUT),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
    CALL timestart("IO (dos-read)")
    SELECT CASE (eig66_data_mode(id))
    CASE (da_mode)
       CALL read_dos_DA(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (hdf_mode)
       CALL read_dos_HDF(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (mem_mode)
       CALL read_dos_Mem(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (MPI_mode)
       CALL read_dos_MPI(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    CASE (-1)
       CALL juDFT_error("Could not read eig-file before opening")
    END SELECT
    CALL timestop("IO (dos-read)")
  END SUBROUTINE read_dos

END MODULE m_eig66_io
