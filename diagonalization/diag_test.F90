!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
PROGRAM diag_test
  USE m_judft
  USE m_types_mat
  USE m_types_mpimat
  USE m_eigen_diag
  USE m_io_matrix
  IMPLICIT NONE

  INTEGER            :: matsize,ne,mode,fid
  INTEGER            :: err,isize
  CHARACTER(len=50)  :: filename
  CLASS(t_mat),ALLOCATABLE :: hmat,smat,ev
  REAL,ALLOCATABLE   :: eig(:)
  LOGICAL            :: l_exist,l_real
  REAL               :: t1,t2
#ifdef CPP_MPI
  Include 'mpif.h'
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,err)
  IF (isize>1) THEN
     ALLOCATE(t_mpimat::hmat)
     ALLOCATE(t_mpimat::smat)
     ALLOCATE(hmat%blacsdata)
     smat%blacsdata=>hmat%blacsdata
     hmat%blacsdata%mpi_comm=MPI_COMM_WORLD
  END IF
#endif
  IF (.NOT.ALLOCATED(hmat)) THEN
     ALLOCATE(t_mat::hmat)
     ALLOCATE(t_mat::smat)
  ENDIF


  ! get filename
  filename=judft_string_for_argument("-file")
  INQUIRE(file=trim(filename)//".hdf",exist=l_exist)
  IF (.NOT.l_exist) CALL judft_error("File specified does not exist")
  

  !l_real,matsize is actually only needed if file is created
  fid=open_matrix(l_real,matsize,2,2,trim(filename))

  CALL read_matrix(hmat,1,fid)
  CALL read_matrix(smat,2,fid)
  SELECT TYPE(hmat)
  TYPE is (t_mpimat)
     SELECT TYPE(smat)
     TYPE is (t_mpimat)
        smat%blacsdata=>hmat%blacsdata!make sure we use same blacs-grids
     END SELECT
     ne=0.15*hmat%global_size1
     ALLOCATE(eig(hmat%global_size1))
  CLASS default
     ne=0.15*hmat%matsize1
     ALLOCATE(eig(hmat%matsize1))
  END SELECT
  
  CALL cpu_TIME(t1)
  mode=0
  CALL eigen_diag(mode,hmat,smat,ne,eig,ev)
  CALL cpu_TIME(t2)
  PRINT *,"No of eigenvalues:",ne
  PRINT *,eig(:ne)

  PRINT *,"Time used:",t1-t2

  CALL close_matrix(fid)

END PROGRAM diag_test
  
  
