!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eig66_hdf
#include "juDFT_env.h"
  !*****************************************************************
  ! DESC:Module for hdf-io of eig-file
  !      To be compatible with f90 interface of HDF, use kind for vars
  !
  !      !ATTENTION before calling openeig and after calling closeeig!
  !      !the hdf library has to be initialized or finalized, respectively
  !
  !      CONTAINS the following subroutines:
  !      openeig        opens file
  !      closeeig       closes file
  !      read_keb       reads kpt, enpara and basis data
  !      read_neig      read no of eigenvalues (and eigenvalues itself)
  !      read_eig       reads eigenvectors
  !      writeeig       saves all data for kpt
  !      writesingleeig saves data for one kpt and energy
  !
  !
  !                          Daniel Wortmann, Tue Nov  512:07:522002
  !*****************************************************************
  USE m_eig66_data
  USE m_types
#ifdef CPP_HDF
  USE hdf5
  USE m_hdf_tools
  IMPLICIT NONE

  PRIVATE
  INTEGER, PARAMETER :: one=1,two=2,three=3,zero=0
  !to have the correct
  !type for array constructors

#endif
  PUBLIC open_eig,close_eig
  PUBLIC read_eig,read_dos,write_dos
  PUBLIC write_eig!,writesingleeig,writeeigc,writebas

CONTAINS
  SUBROUTINE priv_find_data(id,d)
    INTEGER,INTENT(IN)::id
    TYPE(t_data_hdf),POINTER:: d

    CLASS(t_data),POINTER   ::dp
    CALL eig66_find_data(dp,id)
    SELECT TYPE(dp)
    TYPE is (t_data_hdf)
       d=>dp
       CLASS default
       CALL judft_error("BUG: wrong datatype in eig66_hdf")
    END SELECT
  END SUBROUTINE priv_find_data
  !----------------------------------------------------------------------
  SUBROUTINE open_eig(id,mpi_comm,nmat,neig,nkpts,jspins,lmax,nlo,ntype,create,l_real,l_soc,nlotot,readonly,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)

    !*****************************************************************
    !     opens hdf-file for eigenvectors+values
    !*****************************************************************
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: id,mpi_comm
    INTEGER, INTENT(IN) :: nmat,neig,nkpts,jspins,nlo,ntype,lmax,nlotot
    LOGICAL, INTENT(IN) :: create,readonly,l_real,l_soc
    LOGICAL, INTENT(IN),OPTIONAL ::l_dos,l_mcd,l_orb
    CHARACTER(LEN=*),OPTIONAL :: filename
    INTEGER,INTENT(IN),OPTIONAL :: layers,nstars,ncored,nsld,nat

#ifdef CPP_HDF

    INTEGER         :: hdferr,access_mode
    INTEGER(HID_T)  :: creation_prp,access_prp,spaceid
    LOGICAL         :: l_exist
    INTEGER(HSIZE_T):: dims(7)
    TYPE(t_data_HDF),POINTER::d
    !Set creation and access properties
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    IF (readonly) THEN
       access_prp=H5P_DEFAULT_f
       creation_prp=H5P_DEFAULT_f
    ELSE
       CALL h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)
       !      CALL h5pset_fapl_mpiposix_f(access_prp,MPI_COMM,
       !     +.false.,hdferr)
       CALL h5pset_fapl_mpio_f(access_prp, MPI_COMM, MPI_INFO_NULL,hdferr)
       creation_prp=H5P_DEFAULT_f !no special creation property
    ENDIF
#else
    access_prp=H5P_DEFAULT_f
    creation_prp=H5P_DEFAULT_f
#endif 
    CALL priv_find_data(id,d)
    IF (PRESENT(filename)) d%fname=filename
    CALL eig66_data_storedefault(d,jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype,l_real,l_soc,l_dos,l_mcd,l_orb)
    !set access_flags according
    IF (readonly) THEN
       access_mode=H5F_ACC_RDONLY_F
    ELSE
       access_mode=H5F_ACC_RDWR_F
    ENDIF
    !     OPEN FILE and get D%FID's
    IF (create) THEN
       INQUIRE(FILE=TRIM(d%fname)//'.hdf',EXIST=l_exist)
       access_mode=H5F_ACC_TRUNC_F
       !         IF (l_exist) WRITE (*,*)'Warning: eig.hdf was overwritten'
       CALL h5fcreate_f(TRIM(d%fname)//'.hdf',access_Mode, d%fid, hdferr ,creation_prp,access_prp)
       ! create dataspaces and datasets
       !   scalars
       dims(:2)=(/nkpts,jspins/)
       CALL h5screate_simple_f(2,dims(:2),spaceid,hdferr)
       CALL h5dcreate_f(d%fid, "neig", H5T_NATIVE_INTEGER, spaceid, d%neigsetid, hdferr)
       CALL h5dcreate_f(d%fid, "wk", H5T_NATIVE_DOUBLE,spaceid, d%wksetid, hdferr)
       CALL h5dcreate_f(d%fid, "nv", H5T_NATIVE_INTEGER, spaceid, d%nvsetid, hdferr)
       CALL h5dcreate_f(d%fid, "nmat", H5T_NATIVE_INTEGER, spaceid, d%nmatsetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       !   vectors
       dims(1:3)=(/two,nkpts,jspins/)
       CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
       CALL h5dcreate_f(d%fid, "evac", H5T_NATIVE_DOUBLE, spaceid, d%evacsetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       dims(:3)=(/three,nkpts,jspins/)
       CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
       CALL h5dcreate_f(d%fid, "bk", H5T_NATIVE_DOUBLE, spaceid, d%bksetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       dims(:3)=(/neig,nkpts,jspins/)
       CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
       !     ew
       CALL h5dcreate_f(d%fid, "energy", H5T_NATIVE_DOUBLE, spaceid, d%energysetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       !     w_iks
       CALL h5dcreate_f(d%fid, "w_iks", H5T_NATIVE_DOUBLE, spaceid, d%wikssetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       !     enparas
       dims(1:4)=(/lmax+1,ntype,nkpts,jspins/)
       CALL h5screate_simple_f(4,dims(1:4),spaceid,hdferr)
       CALL h5dcreate_f(d%fid, "el", H5T_NATIVE_DOUBLE, spaceid, d%esetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)

       dims(:4)=(/nlo,ntype,nkpts,jspins/)
       CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
       CALL h5dcreate_f(d%fid, "ello", H5T_NATIVE_DOUBLE, spaceid, d%ellosetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       !     ev
       if ( l_real .and..not.l_soc ) THEN
          dims(:5)=(/one,nmat,neig,nkpts,jspins/)
       else
          dims(:5)=(/two,nmat,neig,nkpts,jspins/)
       endif
       CALL h5screate_simple_f(5,dims(:5),spaceid,hdferr)
       CALL h5dcreate_f(d%fid, "ev", H5T_NATIVE_DOUBLE, spaceid, d%evsetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       !      basis
       dims(:4)=(/nmat,three,nkpts,jspins/)
       CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
       CALL h5dcreate_f(d%fid, "k", H5T_NATIVE_INTEGER, spaceid, d%ksetid, hdferr)
       CALL h5sclose_f(spaceid,hdferr)
       !stuff for dos etc
       IF (d%l_dos) THEN
          dims(:5)=(/4,ntype,neig,nkpts,jspins/)
          CALL h5screate_simple_f(5,dims(:5),spaceid,hdferr)
          CALL h5dcreate_f(d%fid, "qal", H5T_NATIVE_DOUBLE, spaceid, d%qalsetid, hdferr)
          CALL h5sclose_f(spaceid,hdferr)
          dims(:4)=(/neig,2,nkpts,jspins/)
          CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
          CALL h5dcreate_f(d%fid, "qvac", H5T_NATIVE_DOUBLE, spaceid, d%qvacsetid, hdferr)
          CALL h5sclose_f(spaceid,hdferr)
          dims(:3)=(/neig,nkpts,jspins/)
          CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
          CALL h5dcreate_f(d%fid, "qis", H5T_NATIVE_DOUBLE, spaceid, d%qissetid, hdferr)
          CALL h5sclose_f(spaceid,hdferr)
          dims(:5)=(/neig,layers,2,nkpts,jspins/)
          CALL h5screate_simple_f(5,dims(:5),spaceid,hdferr)
          CALL h5dcreate_f(d%fid, "qvlay", H5T_NATIVE_DOUBLE, spaceid, d%qvlaysetid, hdferr)
          CALL h5sclose_f(spaceid,hdferr)
          dims(:7)=(/2,nstars,neig,layers,2,nkpts,jspins/)
          CALL h5screate_simple_f(7,dims(:7),spaceid,hdferr)
          CALL h5dcreate_f(d%fid, "qstars", H5T_NATIVE_DOUBLE, spaceid, d%qstarssetid, hdferr)
          CALL h5sclose_f(spaceid,hdferr)
          dims(:3)=(/neig,nkpts,jspins/)
          CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
          CALL h5dcreate_f(d%fid, "ksym", H5T_NATIVE_DOUBLE, spaceid, d%ksymsetid, hdferr)
          CALL h5sclose_f(spaceid,hdferr)
          dims(:3)=(/neig,nkpts,jspins/)
          CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
          CALL h5dcreate_f(d%fid, "jsym", H5T_NATIVE_DOUBLE, spaceid, d%jsymsetid, hdferr)
          CALL h5sclose_f(spaceid,hdferr)
          IF (d%l_mcd) THEN
             dims(:5)=(/3*ntype,ncored,neig,nkpts,jspins/)
             CALL h5screate_simple_f(5,dims(:5),spaceid,hdferr)
             CALL h5dcreate_f(d%fid, "mcd", H5T_NATIVE_DOUBLE, spaceid, d%mcdsetid, hdferr)
             CALL h5sclose_f(spaceid,hdferr)
          ENDIF
          IF (d%l_orb) THEN
             dims(:4)=(/nsld,neig,nkpts,jspins/)
             CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
             CALL h5dcreate_f(d%fid, "qintsl", H5T_NATIVE_DOUBLE, spaceid, d%qintslsetid, hdferr)
             CALL h5sclose_f(spaceid,hdferr)
             dims(:4)=(/nsld,neig,nkpts,jspins/)
             CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
             CALL h5dcreate_f(d%fid, "qmtsl", H5T_NATIVE_DOUBLE, spaceid, d%qmtslsetid, hdferr)
             CALL h5sclose_f(spaceid,hdferr)
             dims(:4)=(/neig,nat,nkpts,jspins/)
             CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
             CALL h5dcreate_f(d%fid, "qmtp", H5T_NATIVE_DOUBLE, spaceid, d%qmtpsetid, hdferr)
             CALL h5sclose_f(spaceid,hdferr)
             dims(:5)=(/neig,23,nat,nkpts,jspins/)
             CALL h5screate_simple_f(5,dims(:5),spaceid,hdferr)
             CALL h5dcreate_f(d%fid, "orbcomp", H5T_NATIVE_DOUBLE, spaceid, d%orbcompsetid, hdferr)
             CALL h5sclose_f(spaceid,hdferr)
          ENDIF
       ENDIF
    ELSE
       CALL h5fopen_f (TRIM(d%fname)//'.hdf', access_Mode, d%fid, hdferr,access_prp)
       !get dataset-ids
       CALL h5dopen_f(d%fid, 'el', d%esetid, hdferr)
       CALL h5dopen_f(d%fid, 'evac', d%evacsetid, hdferr)
       CALL h5dopen_f(d%fid, 'ello', d%ellosetid, hdferr)
       CALL h5dopen_f(d%fid, 'bk', d%bksetid, hdferr)
       CALL h5dopen_f(d%fid, 'wk', d%wksetid, hdferr)
       CALL h5dopen_f(d%fid, 'energy', d%energysetid, hdferr)
       CALL h5dopen_f(d%fid, 'w_iks', d%wikssetid, hdferr)
       CALL h5dopen_f(d%fid, 'k', d%ksetid, hdferr)
       CALL h5dopen_f(d%fid, 'neig', d%neigsetid, hdferr)
       CALL h5dopen_f(d%fid, 'ev', d%evsetid, hdferr)
       CALL h5dopen_f(d%fid, 'nv', d%nvsetid, hdferr)
       CALL h5dopen_f(d%fid, 'nmat', d%nmatsetid, hdferr)
       IF (d%l_dos) THEN
          CALL h5dopen_f(d%fid, 'qal', d%qalsetid, hdferr)
          CALL h5dopen_f(d%fid, 'qvac', d%qvacsetid, hdferr)
          CALL h5dopen_f(d%fid, 'qis', d%qissetid, hdferr)
          CALL h5dopen_f(d%fid, 'qvlay', d%qvlaysetid, hdferr)
          CALL h5dopen_f(d%fid, 'qstars', d%qstarssetid, hdferr)
          CALL h5dopen_f(d%fid, 'ksym', d%ksymsetid, hdferr)
          CALL h5dopen_f(d%fid, 'jsym', d%jsymsetid, hdferr)
          IF (d%l_mcd) THEN
             CALL h5dopen_f(d%fid, 'mcd', d%mcdsetid, hdferr)
          ENDIF
          IF (d%l_orb) THEN
             CALL h5dopen_f(d%fid, 'qintsl', d%qintslsetid, hdferr)
             CALL h5dopen_f(d%fid, 'qmtsl', d%qmtslsetid, hdferr)
             CALL h5dopen_f(d%fid, 'qmtp', d%qmtpsetid, hdferr)
             CALL h5dopen_f(d%fid, 'orbcomp', d%orbcompsetid, hdferr)
          ENDIF
       ENDIF
    endif
    IF (.NOT.access_prp==H5P_DEFAULT_f) CALL H5Pclose_f(access_prp&
            &     ,hdferr)
#else
    CALL juDFT_error("Could not use HDF5 for IO, please recompile")
#endif
  END SUBROUTINE open_eig
     !----------------------------------------------------------------------
  SUBROUTINE close_eig(id,filename)
       !*****************************************************************
       !     closes hdf-file for eigenvectors+values
       !*****************************************************************
       IMPLICIT NONE
       INTEGER,INTENT(IN)                   :: id
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filename

       INTEGER::hdferr
       TYPE(t_data_HDF),POINTER::d

       !close datasets
#ifdef CPP_HDF
       CALL priv_find_data(id,d)

       CALL h5dclose_f(d%esetid,hdferr)
       CALL h5dclose_f(d%evacsetid,hdferr)
       CALL h5dclose_f(d%ellosetid,hdferr)
       CALL h5dclose_f(d%bksetid,hdferr)
       CALL h5dclose_f(d%wksetid,hdferr)
       CALL h5dclose_f(d%energysetid,hdferr)
       CALL h5dclose_f(d%wikssetid,hdferr)
       CALL h5dclose_f(d%ksetid,hdferr)
       CALL h5dclose_f(d%neigsetid,hdferr)
       CALL h5dclose_f(d%evsetid,hdferr)
       CALL h5dclose_f(d%nvsetid,hdferr)
       CALL h5dclose_f(d%nmatsetid,hdferr)
       IF (d%l_dos) THEN
          CALL h5dclose_f(d%qalsetid, hdferr)
          CALL h5dclose_f(d%qvacsetid, hdferr)
          CALL h5dclose_f(d%qissetid, hdferr)
          CALL h5dclose_f(d%qvlaysetid, hdferr)
          CALL h5dclose_f(d%qstarssetid, hdferr)
          CALL h5dclose_f(d%ksymsetid, hdferr)
          CALL h5dclose_f(d%jsymsetid, hdferr)
          IF (d%l_mcd) THEN
             CALL h5dclose_f(d%mcdsetid, hdferr)
          ENDIF
          IF (d%l_orb) THEN
             CALL h5dclose_f(d%qintslsetid, hdferr)
             CALL h5dclose_f(d%qmtslsetid, hdferr)
             CALL h5dclose_f(d%qmtpsetid, hdferr)
             CALL h5dclose_f(d%orbcompsetid, hdferr)
          ENDIF
       ENDIF
       !close file
       CALL h5fclose_f(d%fid,hdferr)
       !If a filename was given and the name is not the current filename
       IF (PRESENT(filename)) THEN
          IF (filename.NE.d%fname) THEN
             CALL system("mv "//TRIM(d%fname)//".hdf "//TRIM(filename)//".hdf")
          ENDIF
       ENDIF
       d%fname="eig"
       CALL eig66_remove_data(id)

#endif
     END SUBROUTINE close_eig
#ifdef CPP_HDF
     !----------------------------------------------------------------------
     SUBROUTINE priv_r_vec(d,nk,jspin,n_start,n_end,nmat,z)

       USE m_hdf_tools
       IMPLICIT NONE
       TYPE(t_data_HDF),INTENT(IN)::d
       INTEGER, INTENT(IN)  :: nk,jspin
       INTEGER, INTENT(IN)  :: n_start,n_end
       INTEGER, INTENT(OUT) :: nmat
       REAL,    INTENT(OUT) :: z(:,:)

       INTEGER i,j,neig_l

       neig_l = n_end - n_start + 1

       ! read matrix size
       CALL io_read_integer0(d%nmatsetid,(/nk,jspin/),(/1,1/),nmat)

       IF ( nmat > SIZE(z,1) .OR. neig_l > SIZE(z,2) ) THEN
          WRITE (6,*) nmat,SIZE(z,1),SIZE(z,2)
          CALL juDFT_error("eig66_hdf$read_vec",calledby ="eig66_hdf")
       ENDIF

       !read eigenvectors
       CALL io_read_real2(d%evsetid,(/1,1,n_start,nk,jspin/),&
            &                           (/1,nmat,neig_l,1,1/),&
            &                           z(:nmat,:neig_l) )

     END SUBROUTINE priv_r_vec

#endif
     SUBROUTINE read_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
       IMPLICIT NONE
       INTEGER, INTENT(IN)          :: id,nk,jspin
       REAL,INTENT(OUT)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
       COMPLEX,INTENT(OUT)           :: qstars(:,:,:,:)
       INTEGER,INTENT(OUT)           :: ksym(:),jsym(:)
       REAL,INTENT(OUT),OPTIONAL     :: mcd(:,:,:)
       REAL,INTENT(OUT),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
       TYPE(t_data_HDF),POINTER      :: d
       REAL,ALLOCATABLE              :: r_tmp5(:,:,:,:,:)
       CALL priv_find_data(id,d)
#ifdef CPP_HDF
       CALL io_read_real3(d%qalsetid,(/1,1,1,nk,jspin/),(/SIZE(qal,1),SIZE(qal,2),SIZE(qal,3),1,1/),qal)
       CALL io_read_real2(d%qvacsetid,(/1,1,nk,jspin/),(/SIZE(qvac,1),SIZE(qvac,2),1,1/),qvac)
       CALL io_read_real1(d%qissetid,(/1,nk,jspin/),(/SIZE(qis,1),1,1/),qis)
       CALL io_read_real3(d%qvlaysetid,(/1,1,1,nk,jspin/),(/SIZE(qvlay,1),SIZE(qvlay,2),SIZE(qvlay,3),1,1/),qvlay)
       ALLOCATE(r_tmp5(2,SIZE(qstars,1),SIZE(qstars,2),SIZE(qstars,3),SIZE(qstars,4)))
       CALL io_read_real5(d%qstarssetid,(/1,1,1,1,1,nk,jspin/),(/2,SIZE(qstars,1),SIZE(qstars,2),SIZE(qstars,3),SIZE(qstars,4),1,1/),r_tmp5(:,:,:,:,:))
       qstars=CMPLX(r_tmp5(1,:,:,:,:),r_tmp5(2,:,:,:,:))
       DEALLOCATE(r_tmp5)
       CALL io_read_integer1(d%ksymsetid,(/1,nk,jspin/),(/SIZE(ksym,1),1,1/),ksym)
       CALL io_read_integer1(d%jsymsetid,(/1,nk,jspin/),(/SIZE(jsym,1),1,1/),jsym)
       IF (d%l_mcd.AND.PRESENT(mcd)) THEN
          CALL io_read_real3(d%mcdsetid,(/1,1,1,nk,jspin/),(/SIZE(mcd,1),SIZE(mcd,2),SIZE(mcd,3),1,1/),mcd)
       ENDIF
       IF (d%l_orb.AND.PRESENT(qintsl)) THEN
          CALL io_read_real2(d%qintslsetid,(/1,1,nk,jspin/),(/SIZE(qintsl,1),SIZE(qintsl,2),1,1/),qintsl)
          CALL io_read_real2(d%qmtslsetid,(/1,1,nk,jspin/),(/SIZE(qmtsl,1),SIZE(qmtsl,2),1,1/),qmtsl)
          CALL io_read_real2(d%qmtpsetid,(/1,1,nk,jspin/),(/SIZE(qmtp,1),SIZE(qmtp,2),1,1/),qmtp)
          CALL io_read_real3(d%orbcompsetid,(/1,1,1,nk,jspin/),(/SIZE(orbcomp,1),23,SIZE(orbcomp,3),1,1/),orbcomp)
       ENDIF
#endif
     END SUBROUTINE read_dos


     SUBROUTINE write_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
       IMPLICIT NONE
       INTEGER, INTENT(IN)          :: id,nk,jspin
       REAL,INTENT(IN)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
       COMPLEX,INTENT(IN)           :: qstars(:,:,:,:)
       INTEGER,INTENT(IN)           :: ksym(:),jsym(:)
       REAL,INTENT(IN),OPTIONAL     :: mcd(:,:,:)
       REAL,INTENT(IN),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
       TYPE(t_data_HDF),POINTER      ::d
       CALL priv_find_data(id,d)
#ifdef CPP_HDF
       CALL io_write_real3(d%qalsetid,(/1,1,1,nk,jspin/),(/SIZE(qal,1),SIZE(qal,2),SIZE(qal,3),1,1/),qal)
       CALL io_write_real2(d%qvacsetid,(/1,1,nk,jspin/),(/SIZE(qvac,1),SIZE(qvac,2),1,1/),qvac)
       CALL io_write_real1(d%qissetid,(/1,nk,jspin/),(/SIZE(qis,1),1,1/),qis)
       CALL io_write_real3(d%qvlaysetid,(/1,1,1,nk,jspin/),(/SIZE(qvlay,1),SIZE(qvlay,2),SIZE(qvlay,3),1,1/),qvlay)
       CALL io_write_real4(d%qstarssetid,(/1,1,1,1,1,nk,jspin/),(/1,SIZE(qstars,1),SIZE(qstars,2),SIZE(qstars,3),SIZE(qstars,4),1,1/),REAL(qstars))
       CALL io_write_real4(d%qstarssetid,(/2,1,1,1,1,nk,jspin/),(/1,SIZE(qstars,1),SIZE(qstars,2),SIZE(qstars,3),SIZE(qstars,4),1,1/),AIMAG(qstars))

       CALL io_write_integer1(d%ksymsetid,(/1,nk,jspin/),(/SIZE(ksym,1),1,1/),ksym)
       CALL io_write_integer1(d%jsymsetid,(/1,nk,jspin/),(/SIZE(jsym,1),1,1/),jsym)
       IF (d%l_mcd.AND.PRESENT(mcd)) THEN
          CALL io_write_real3(d%mcdsetid,(/1,1,1,nk,jspin/),(/SIZE(mcd,1),SIZE(mcd,2),SIZE(mcd,3),1,1/),mcd)
       ENDIF
       IF (d%l_orb.AND.PRESENT(qintsl)) THEN
          CALL io_write_real2(d%qintslsetid,(/1,1,nk,jspin/),(/SIZE(qintsl,1),SIZE(qintsl,2),1,1/),qintsl)
          CALL io_write_real2(d%qmtslsetid,(/1,1,nk,jspin/),(/SIZE(qmtsl,1),SIZE(qmtsl,2),1,1/),qmtsl)
          CALL io_write_real2(d%qmtpsetid,(/1,1,nk,jspin/),(/SIZE(qmtp,1),SIZE(qmtp,2),1,1/),qmtp)
          CALL io_write_real3(d%orbcompsetid,(/1,1,1,nk,jspin/),(/SIZE(orbcomp,1),23,SIZE(orbcomp,3),1,1/),orbcomp)
       ENDIF
#endif
     END SUBROUTINE write_dos


     SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,bk,wk,&
          &                  eig,w_iks,el,ello,evac,&
          &                  nlotot,n_size,n_rank,zmat)

       !*****************************************************************
       !     writes all eignevecs for the nk-th kpoint
       !*****************************************************************
       IMPLICIT NONE

       INTEGER, INTENT(IN)          :: id,nk,jspin
       INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
       REAL,    INTENT(IN),OPTIONAL :: wk
       INTEGER, INTENT(IN),OPTIONAL :: neig,nv,nmat,nlotot,neig_total
       REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:),w_iks(:)
       REAL,    INTENT(IN),OPTIONAL :: evac(2),ello(:,:)
       TYPE(t_zmat),INTENT(IN),OPTIONAL :: zmat

       INTEGER i,j,k,nv_local,n1,n2,ne
       TYPE(t_data_HDF),POINTER::d
       CALL priv_find_data(id,d)

#ifdef CPP_HDF
       !
       !write enparas
       !
       nv_local=HUGE(1)

       IF (PRESENT(el))&
            &   CALL io_write_real2(&
            &                    d%esetid,(/1,1,nk,jspin/),&
            &                    (/SIZE(el,1),SIZE(el,2),1,1/),el)

       IF (PRESENT(ello))&
            & CALL io_write_real2(&
            &                    d%ellosetid,(/1,1,nk,jspin/),&
            &                    (/SIZE(ello,1),SIZE(ello,2),1,1/),ello)

       IF (PRESENT(evac)) CALL io_write_real1(&
            &                    d%evacsetid,(/1,nk,jspin/),(/2,1,1/),evac)
       !
       !write kpts
       !

       IF (PRESENT(bk)) CALL io_write_real1(&
            &                    d%bksetid,(/1,nk,jspin/),(/3,1,1/),bk)

       IF (PRESENT(wk)) CALL io_write_real0(&
            &                    d%wksetid,(/nk,jspin/),(/1,1/),wk)
       !
       !write basis
       !

       IF (PRESENT(nv)) THEN
          nv_local=nv
          CALL io_write_integer0(d%nvsetid,(/nk,jspin/),(/1,1/),nv)
       ENDIF

       IF (PRESENT(nmat)) CALL io_write_integer0(&
            &                       d%nmatsetid,(/nk,jspin/),(/1,1/),nmat)

       ENDIF
       !
       !write eigenvalues
       !
       IF (PRESENT(w_iks)) THEN
          CALL io_write_real1s(d%wikssetid,(/1,nk,jspin/),(/size(w_iks),1,1/),w_iks,(/1,1,1/))
       ENDIF
       
       IF (PRESENT(neig_total)) THEN
          CALL io_write_integer0(d%neigsetid,(/nk,jspin/),(/1,1/),neig_total)
       ENDIF

       IF (PRESENT(n_rank).AND.PRESENT(n_size).AND.&
            &        PRESENT(eig).AND.PRESENT(neig)) THEN
          CALL io_write_real1s(&
               &                     d%energysetid,(/n_rank+1,nk,jspin/),        &
               &                     (/neig,1,1/),eig(:neig),(/n_size,1,1/))
          !write eigenvectors
          !
       ELSEIF (PRESENT(eig).AND.PRESENT(neig)) THEN
          CALL io_write_real1s(&
               &                     d%energysetid,(/1,nk,jspin/),&
               &                     (/neig,1,1/),eig(:neig),(/1,1,1/))
       ELSE
          IF (PRESENT(eig)) CALL juDFT_error("BUG in calling write_eig")
       ENDIF
       IF (PRESENT(zmat).AND..NOT.PRESENT(neig))&
            &    CALL juDFT_error("BUG in calling write_eig with eigenvector")

       n1=1;n2=0
       IF (PRESENT(n_size)) n1=n_size
       IF (PRESENT(n_rank)) n2=n_rank
       IF (PRESENT(zmat)) THEN
          IF (zmat%l_real) THEN
             CALL io_write_real2s(&
                  &                     d%evsetid,(/1,1,n2+1,nk,jspin/),&
                  &           (/1,nmat,neig,1,1/),REAL(zmat%z_r(:nmat,:neig)),(/1,1,n1,1,1/))
          ELSE
             CALL io_write_real2s(&
                  &                     d%evsetid,(/1,1,n2+1,nk,jspin/),&
                  &           (/1,nmat,neig,1,1/),REAL(zmat%z_c(:nmat,:neig)),(/1,1,n1,1,1/))
             CALL io_write_real2s(&
                  &                     d%evsetid,(/2,1,n2+1,nk,jspin/),&
                  &           (/1,nmat,neig,1,1/),AIMAG(zmat%z_c(:nmat,:neig)),&
                  &           (/1,1,n1,1,1/))
          ENDIF
       ENDIF

#endif
     END SUBROUTINE write_eig

#ifdef CPP_HDF

     !----------------------------------------------------------------------
     SUBROUTINE priv_r_vecc(&
          &                     d,nk,jspin,n_start,n_end,&
          &                     nmat,z)

       USE m_hdf_tools
       IMPLICIT NONE
       TYPE(t_data_HDF),INTENT(IN)::d
       INTEGER, INTENT(IN)  :: nk,jspin
       INTEGER, INTENT(IN)  :: n_start,n_end
       INTEGER, INTENT(OUT) :: nmat
       COMPLEX, INTENT(OUT) :: z(:,:)

       REAL, ALLOCATABLE :: z1(:,:,:)
       INTEGER i,j,neig_l

       neig_l = n_end - n_start + 1

       ! read matrix size
       CALL io_read_integer0(&
            &                      d%nmatsetid,(/nk,jspin/),(/1,1/),&
            &                                                nmat)

       IF ( nmat > SIZE(z,1) .OR. neig_l > SIZE(z,2) ) THEN
          WRITE (6,*) nmat,SIZE(z,1),SIZE(z,2)
          CALL juDFT_error("eig66_hdf$read_vec",calledby ="eig66_hdf")
       ENDIF

       ! read eigenvectors
       ALLOCATE (z1(2,nmat,neig_l))
       CALL io_read_real3(d%evsetid,(/1,1,n_start,nk,jspin/),&
            &                      (/2,nmat,neig_l,1,1/),z1)

       DO i=1,neig_l
          DO j=1,nmat
             z(j,i) = CMPLX( z1(1,j,i) ,z1(2,j,i) )
          ENDDO
       ENDDO

       DEALLOCATE (z1)

     END SUBROUTINE priv_r_vecc
     !-----------------------------------------------------------------------

#endif

     SUBROUTINE read_eig(id,nk,jspin,nv,nmat,bk,wk,neig,eig,w_iks,el,&
          &            ello,evac,n_start,n_end,zMat)
       IMPLICIT NONE
       INTEGER, INTENT(IN)            :: id,nk,jspin
       INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
       INTEGER, INTENT(OUT),OPTIONAL  :: neig
       REAL,    INTENT(OUT),OPTIONAL  :: eig(:),w_iks(:)
       REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
       REAL,    INTENT(OUT),OPTIONAL  :: bk(:),wk
       INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
       TYPE(t_zMat),OPTIONAL  :: zmat

#ifdef CPP_HDF
       INTEGER:: n1,n,k
       TYPE(t_data_HDF),POINTER::d
       CALL priv_find_data(id,d)


       IF (PRESENT(neig))  THEN
          CALL io_read_integer0(d%neigsetid,(/nk,jspin/),(/1,1/),neig)

          IF ( PRESENT(eig) ) THEN                           ! read eigenv
             IF ( neig > SIZE(eig) ) THEN
                WRITE(*,*) neig,SIZE(eig)
                CALL juDFT_error("eig66_hdf$readeig",calledby ="eig66_hdf")
             ENDIF
             CALL io_read_real1(d%energysetid,(/1,nk,jspin/),(/neig,1,1/),&
                  &                      eig(:neig))
          ENDIF
          IF (PRESENT(w_iks)) THEN
             CALL io_read_real1(d%wikssetid,(/1,nk,jspin/),(/size(w_iks),1,1/),w_iks)
          ENDIF
       ENDIF

     
          IF (PRESENT(nv)) nv=n1
       ELSE
          IF (PRESENT(nv)) CALL io_read_integer0(d%nvsetid,(/nk,jspin/),(/1,1/),nv)

       ENDIF
       IF (PRESENT(nmat)) &
            & CALL io_read_integer0(d%nmatsetid,(/nk,jspin/),(/1,1/),nmat)
       IF (PRESENT(el)) CALL io_read_real2(d%esetid,(/1,1,nk,jspin/),&
            &                   (/SIZE(el,1),SIZE(el,2),1,1/),el(:,:))
       IF (PRESENT(ello)) CALL io_read_real2(d%ellosetid,(/1,1,nk,jspin/),&
            &                   (/SIZE(ello,1),SIZE(ello,2),1,1/),ello(:,:))
       IF (PRESENT(evac)) CALL io_read_real1(d%evacsetid,(/1,nk,jspin/),&
            &                 (/2,1,1/),evac)

       IF (PRESENT(bk)) CALL&
            io_read_real1(d%bksetid,(/1,nk,jspin/),(/3,1,1/),bk)
       IF (PRESENT(wk)) CALL&
            io_read_real0(d%wksetid,(/nk,jspin/),(/1,1/),wk)

       IF (PRESENT(n_start)) THEN
          IF (.NOT.PRESENT(n_end)) CALL juDFT_error("BUG3 in read_eig")
          IF (PRESENT(zMat)) THEN
             IF (zmat%l_real) THEN
                CALL priv_r_vec(d,nk,jspin,n_start,n_end,n1,zmat%z_r)
             ELSE
                CALL priv_r_vecc(d,nk,jspin,n_start,n_end,n1,zmat%z_c)
             ENDIF
          ENDIF
          IF (PRESENT(nmat)) nmat=n1
       ENDIF
#endif
     END SUBROUTINE read_eig

   END MODULE

