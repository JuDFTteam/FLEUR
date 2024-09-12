!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eig66_io
use m_juDFT
   use m_types_mat
   USE m_eig66_data
   IMPLICIT NONE
   PRIVATE

   PUBLIC open_eig, close_eig, reset_eig
   PUBLIC read_eig, write_eig
CONTAINS

   FUNCTION open_eig(mpi_comm, nmat, neig, nkpts, jspins, &
                     l_noco, l_create, l_real, l_soc, l_readonly, l_olap, n_size, mode_in, filename) &
      RESULT(id)
      USE m_eig66_hdf, ONLY: open_eig_hdf => open_eig
      USE m_eig66_DA, ONLY: open_eig_DA => open_eig
      USE m_eig66_mem, ONLY: open_eig_mem => open_eig
      USE m_eig66_MPI, ONLY: open_eig_mpi => open_eig
      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: nmat, neig, nkpts, jspins, mpi_comm
      LOGICAL, INTENT(IN)          :: l_noco, l_readonly, l_create, l_real, l_soc, l_olap
      INTEGER, INTENT(IN), OPTIONAL :: n_size, mode_in
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
      INTEGER:: id, mode
      CHARACTER(LEN=20) :: arg

      INTEGER:: neig_local, isize, err
      if (l_soc) THEN
         neig_local = 2*neig
      else
         neig_local = neig
      endif
      mode = -1
      IF (PRESENT(mode_in)) mode = mode_in

      IF (mode < 0) THEN
         !Use default mode
#ifdef CPP_MPI
         CALL MPI_COMM_SIZE(mpi_comm, isize, err)
         IF (isize > 1) THEN
            mode = MPI_mode
         ELSE
            mode = MEM_mode
         ENDIF
#else
         mode = MEM_mode
#endif
         !check if default was given on command-line
         arg=TRIM(juDFT_string_for_argument("-eig"))

         IF (index(arg,"mpi")>0) mode = MPI_mode
         IF (index(arg,"mem")>0) mode = MEM_mode
         IF (index(arg,"da")>0) mode = DA_mode
         IF (index(arg,"hdf")>0) mode = HDF_mode
      ENDIF
      !Check if mode is available
#ifndef CPP_MPI
      IF (mode == MPI_mode) CALL juDFT_error("MPI-mode not available. Recompile with CPP_MPI", calledby="eig66_io")
#else
      CALL MPI_COMM_SIZE(mpi_comm, isize, err)
      IF (isize > 1 .AND. ((mode == DA_mode .OR. mode == mem_mode))) &
         CALL juDFT_error("In a parallel calculation MEM/DA-mode are not available", calledby="eig66_io")
#endif
#ifndef CPP_HDF
      IF (mode == HDF_mode) CALL juDFT_error("HDF-mode not available. Recompile with CPP_HDF", calledby="eig66_io")
#endif

      id = eig66_data_newid(mode)

      !PRINT *,"open_eig:",id,mode

      CALL timestart("Open file/memory for IO of eig66")
      SELECT CASE (eig66_data_mode(id))
      CASE (DA_mode)
         CALL open_eig_DA(id, nmat, neig_local, nkpts, jspins, l_create, l_real, l_soc, l_olap, filename)
      CASE (hdf_mode)
         CALL open_eig_HDF(id, mpi_comm, nmat, neig_local, nkpts, jspins, l_create, l_real, l_soc, l_readonly, l_olap, filename)
      CASE (mem_mode)
         CALL open_eig_MEM(id, nmat, neig_local, nkpts, jspins, l_create, l_real, l_soc, l_noco, l_olap, filename)
      CASE (mpi_mode)
         CALL open_eig_MPI(id, mpi_comm, nmat, neig_local, nkpts, jspins, l_create, l_real, l_soc, l_noco, l_olap, n_size, filename)
      CASE DEFAULT
         CALL juDFT_error("Invalid IO-mode in eig66_io")
      END SELECT
      CALL timestop("Open file/memory for IO of eig66")
   END FUNCTION open_eig

   SUBROUTINE close_eig(id, filename)
      USE m_eig66_hdf, ONLY: close_eig_hdf => close_eig
      USE m_eig66_DA, ONLY: close_eig_DA => close_eig
      USE m_eig66_mem, ONLY: close_eig_MEM => close_eig
      USE m_eig66_MPI, ONLY: close_eig_MPI => close_eig
      IMPLICIT NONE
      INTEGER, INTENT(IN)                   :: id
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
      INTEGER  :: mode
      mode = eig66_data_mode(id)
      !PRINT*,"close_eig:",id,mode
      SELECT CASE (mode)
      CASE (DA_mode)
         CALL close_eig_DA(id, filename)
      CASE (hdf_mode)
         CALL close_eig_HDF(id, filename)
      CASE (mem_mode)
         CALL close_eig_Mem(id, filename=filename)
      CASE (MPI_mode)
         CALL close_eig_MPI(id, filename=filename)
      CASE (-1)
         CALL juDFT_error("ID not assigned in close_eig", calledby="eig66_io")
      END SELECT

   END SUBROUTINE close_eig

   subroutine read_eig(id, nk, jspin, neig, eig, list, zmat, smat,kpts,input,noco,nococonv,sym,atoms,cell)
      use m_types
      use m_types_nococonv
      use m_trafo
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: id, nk, jspin
      INTEGER, INTENT(OUT), OPTIONAL  :: neig
      REAL, INTENT(OUT), OPTIONAL  :: eig(:)
      INTEGER, INTENT(IN), OPTIONAL   :: list(:)
      TYPE(t_mat), INTENT(INOUT), OPTIONAL  :: zmat, smat
      !for the rotation to IBZ we need:
      TYPE(t_kpts),INTENT(IN),OPTIONAL  :: kpts
      TYPE(t_input),INTENT(IN),OPTIONAL :: input 
      TYPE(t_noco),INTENT(IN),OPTIONAL  :: noco
      TYPE(t_nococonv),INTENT(IN),OPTIONAL:: nococonv
      TYPE(t_sym),INTENT(IN),OPTIONAL   :: sym
      TYPE(t_atoms),INTENT(IN),OPTIONAL :: atoms
      TYPE(t_cell),INTENT(IN),OPTIONAL  :: cell
      

      LOGICAL:: ibz
      TYPE(t_mat):: tmp_mat
      INTEGER :: ikp,iop
      TYPE(t_lapw):: lapw_nk,lapw_ikp
      ibz=.true.
      if (present(kpts)) THEN
         if (nk>kpts%nkpt) ibz=.false.
      endif   


      if(ibz) then
         call read_eig_ibz(id,nk,jspin,neig, eig=eig, list=list,zmat=zmat,smat=smat) 
      else
        if (present(smat)) call judft_error("SMAT cannot be rotated")
        call tmp_mat%init(zmat)
     
        ikp = kpts%bkp(nk) ! parent k-point
        iop = kpts%bksym(nk) ! connecting symm

        call read_eig_ibz(id,ikp, jspin,zmat=tmp_mat, list=list, eig=eig) 
      
        CALL lapw_nk%init(input, noco, nococonv, kpts, atoms, sym, nk, cell)
        CALL lapw_ikp%init(input, noco, nococonv, kpts, atoms, sym, ikp, cell)
        call waveftrafo_gen_zmat(tmp_mat, ikp, iop, kpts, sym, jspin, tmp_mat%matsize2, lapw_ikp, lapw_nk, zmat)
      endif
   end subroutine read_eig

   SUBROUTINE read_eig_ibz(id, nk, jspin, neig, eig, list, zmat, smat)
      USE m_eig66_hdf, ONLY: read_eig_hdf => read_eig
      USE m_eig66_DA, ONLY: read_eig_DA => read_eig
      USE m_eig66_mem, ONLY: read_eig_mem => read_eig
      USE m_eig66_MPI, ONLY: read_eig_MPI => read_eig
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: id, nk, jspin
      INTEGER, INTENT(OUT), OPTIONAL  :: neig
      REAL, INTENT(OUT), OPTIONAL  :: eig(:)
      INTEGER, INTENT(IN), OPTIONAL   :: list(:)
      TYPE(t_mat), INTENT(INOUT), OPTIONAL  :: zmat, smat
      INTEGER::n
      CALL timestart("IO (read)")
      SELECT CASE (eig66_data_mode(id))
      CASE (DA_mode)
         CALL read_eig_DA(id, nk, jspin, neig, eig, list, zmat, smat)
      CASE (hdf_mode)
         CALL read_eig_hdf(id, nk, jspin, neig, eig, list, zmat, smat)
      CASE (mem_mode)
         CALL read_eig_mem(id, nk, jspin, neig, eig, list, zmat, smat)
      CASE (mpi_mode)
         CALL read_eig_mpi(id, nk, jspin, neig, eig, list, zmat, smat)
      CASE (-1)
         CALL juDFT_error("Could not read eig-file before opening", calledby="eig66_io")
      END SELECT
      CALL timestop("IO (read)")
   END SUBROUTINE read_eig_ibz

   SUBROUTINE write_eig(id, nk, jspin, neig, neig_total, eig, n_start, n_end, zmat, smat)
      USE m_eig66_hdf, ONLY: write_eig_hdf => write_eig
      USE m_eig66_DA, ONLY: write_eig_DA => write_eig
      USE m_eig66_mem, ONLY: write_eig_MEM => write_eig
      USE m_eig66_MPI, ONLY: write_eig_MPI => write_eig
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: id, nk, jspin
      INTEGER, INTENT(IN), OPTIONAL :: neig, neig_total, n_start, n_end
      REAL, INTENT(IN), OPTIONAL :: eig(:)
      TYPE(t_Mat), INTENT(IN), OPTIONAL :: zmat, smat
      CALL timestart("IO (write)")
      SELECT CASE (eig66_data_mode(id))
      CASE (da_mode)
         CALL write_eig_DA(id, nk, jspin, neig, neig_total, eig, n_start, n_end, zmat, smat)
      CASE (hdf_mode)
         CALL write_eig_HDF(id, nk, jspin, neig, neig_total, eig, n_start, n_end, zmat, smat)
      CASE (mem_mode)
         CALL write_eig_Mem(id, nk, jspin, neig, neig_total, eig, n_start, n_end, zmat, smat)
      CASE (MPI_mode)
         CALL write_eig_MPI(id, nk, jspin, neig, neig_total, eig, n_start, n_end, zmat, smat)
      CASE (-1)
         CALL juDFT_error("Could not write eig-file before opening", calledby="eig66_io")
      END SELECT
      CALL timestop("IO (write)")
   END SUBROUTINE write_eig

   SUBROUTINE reset_eig(id, l_soc)
      USE m_eig66_MPI, ONLY: reset_eig_MPI => reset_eig
      INTEGER, INTENT(IN) :: id
      LOGICAL, INTENT(IN) :: l_soc

      SELECT CASE (eig66_data_mode(id))
      CASE (MPI_mode)
         CALL reset_eig_MPI(id, l_soc)
      END SELECT

   END SUBROUTINE reset_eig

END MODULE m_eig66_io
