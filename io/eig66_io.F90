MODULE m_eig66_io
#include "juDFT_env.h"
    use m_eig66_data
    IMPLICIT NONE
    PRIVATE

    PUBLIC open_eig,close_eig
    PUBLIC read_eig, write_eig
    CONTAINS

    FUNCTION open_eig(mpi_comm,nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,l_noco,l_create,l_readonly,n_size,mode_in,filename)result(id)
    USE m_eig66_hdf,ONLY:open_eig_hdf=>open_eig
    USE m_eig66_DA ,ONLY:open_eig_DA=>open_eig
    USE m_eig66_mem,ONLY:open_eig_mem=>open_eig
    USE m_eig66_MPI,ONLY:open_eig_mpi=>open_eig
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,mpi_comm
    LOGICAL,INTENT(IN)          :: l_noco,l_readonly,l_create
    INTEGER,INTENT(IN),OPTIONAL :: n_size,mode_in
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
    INTEGER:: id,mode

    INTEGER:: neig_local,isize,err
#ifdef CPP_SOC
    neig_local=2*neig
#else
    neig_local=neig
#endif
    mode=-1
    if (present(mode_in)) mode=mode_in

    IF (mode<0) THEN
       !Use default mode
#ifdef CPP_MPI
       mode=MPI_mode
#else
       mode=MEM_mode
#endif
       !check if default was given on command-line
       if (juDFT_was_argument("-mpi")) mode=MPI_mode
       if (juDFT_was_argument("-mem")) mode=MEM_mode
       if (juDFT_was_argument("-da")) mode=DA_mode
       if (juDFT_was_argument("-hdf")) mode=HDF_mode
    ENDIF
    !Check if mode is available
#ifndef CPP_MPI
     if (mode==MPI_mode) call juDFT_error("MPI-mode not available. Recompile with CPP_MPI",calledby="eig66_io")
#else
     CALL MPI_COMM_SIZE(mpi_comm,isize,err)
     if (isize>1.and.((mode==DA_mode.or.mode==mem_mode))) &
       call juDFT_error("In a parallel calculation MEM/DA-mode are not available",calledby="eig66_io")
#endif
#ifndef CPP_HDF
     if (mode==HDF_mode) call juDFT_error("HDF-mode not available. Recompile with CPP_HDF",calledby="eig66_io")
#endif

     id=eig66_data_newid(mode)

     print *,"open_eig:",id,mode

    call timestart("Open file/memory for IO of eig66")
    SELECT CASE (eig66_data_mode(id))
       CASE (DA_mode)
            CALL open_eig_DA(id,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,nlotot,l_create,filename)
       CASE (hdf_mode)
            CALL open_eig_HDF(id,mpi_comm,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,l_create,l_readonly,filename)
       CASE (mem_mode)
            CALL open_eig_MEM(id,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,l_create,nlotot,l_noco,filename)
       CASE (mpi_mode)
            CALL open_eig_MPI(id,mpi_comm,nmat,neig_local,nkpts,jspins,lmax,nlo,ntype,l_create,nlotot,l_noco,n_size,filename)
       CASE DEFAULT
            CALL juDFT_error("Invalid IO-mode in eig66_io")
    END SELECT
    call timestop("Open file/memory for IO of eig66")
    END FUNCTION

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
    print*,"close_eig:",id,mode
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

    END SUBROUTINE

    SUBROUTINE read_eig(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,z)
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
    CLASS(*),INTENT(OUT),OPTIONAL  :: z(:,:)
    INTEGER::n
    CALL timestart("IO (read)")
    SELECT CASE (eig66_data_mode(id))
       CASE (DA_mode)
            CALL read_eig_DA(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,z)
       CASE (hdf_mode)
            CALL read_eig_hdf(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,z)
       CASE (mem_mode)
            CALL read_eig_mem(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,z)
       CASE (mpi_mode)
            CALL read_eig_mpi(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,z)
       CASE (-1)
            CALL juDFT_error("Could not read eig-file before opening")
    END SELECT
    call timestop("IO (read)")
    END SUBROUTINE

    SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,z)
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
    CLASS(*),INTENT(IN),OPTIONAL :: z(:,:)
    call timestart("IO (write)")
    SELECT CASE (eig66_data_mode(id))
       CASE (da_mode)
            CALL write_eig_DA(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,z)
       CASE (hdf_mode)
            CALL write_eig_HDF(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,z)
       CASE (mem_mode)
            CALL write_eig_Mem(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,z)
       CASE (MPI_mode)
            CALL write_eig_MPI(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,eig,el,ello,evac,nlotot,kveclo,n_start,n_end,z)
       CASE (-1)
            CALL juDFT_error("Could not write eig-file before opening")
    END SELECT
    call timestop("IO (write)")
    END SUBROUTINE

END MODULE m_eig66_io
