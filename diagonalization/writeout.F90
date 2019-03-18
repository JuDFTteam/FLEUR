MODULE m_writeout
CONTAINS

  SUBROUTINE diag_writeout(smat,hmat)
    USE m_types
    USE m_judft
    USE m_io_matrix
    USE m_types_mpimat
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT) :: hmat,smat
    !small subroutine that does only wite the matrix to a file
    INTEGER:: i,ii,irank,ierr,matsize
    CHARACTER(len=20)::filename
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
#else
    irank=0
#endif

    SELECT TYPE(hmat)
    TYPE is (t_mpimat)
       matsize=hmat%global_size1
    CLASS default
       matsize=hmat%matsize1
    END SELECT
    !First write binary file
#ifdef CPP_HDF
    i=open_matrix(hmat%l_real,matsize,2,2,"hs_mat")
#else
    i=open_matrix(hmat%l_real,hmat%matsize1,1,2,"hs_mat")
#endif
    CALL write_matrix(hmat,1,i)
    CALL write_matrix(smat,2,i)
    CALL close_matrix(i)

    !Now the formatted matrix
    WRITE(filename,"(a,i0)") "hmat",irank
    OPEN(999,file=TRIM(filename))
    WRITE(filename,"(a,i0)") "smat",irank
    OPEN(998,file=TRIM(filename))
    DO i=1,hmat%matsize2
       DO ii=1,hmat%matsize1
          IF (hmat%l_real) THEN
             WRITE(999,"(2i6,f15.6)") ii,i,hmat%data_r(ii,i)
             WRITE(998,"(2i6,f15.6)") ii,i,smat%data_r(ii,i)
          ELSE
             WRITE(999,"(2i6,2f15.6)") ii,i,hmat%data_c(ii,i)
             WRITE(998,"(2i6,2f15.6)") ii,i,smat%data_c(ii,i)
          ENDIF
       END DO
    ENDDO
    CLOSE(999)
    CLOSE(998)
#ifdef CPP_MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
    CALL judft_error("STOP in eigen_diag:debug_diag")
  END SUBROUTINE diag_writeout
END MODULE m_writeout
