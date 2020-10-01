!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mingeselle
  USE m_types_mat
  USE m_types_mpimat
  USE m_juDFT
  IMPLICIT NONE
CONTAINS

  SUBROUTINE mingeselle(mat1,mat)
    !---------------------------------------------------------------------+
    !                                                                     |
    ! Transfers the spin-down/spin-up part , upper triangle of the        |
    ! mat1 to the lower triangle of the mat.                              |
    !                                                                     |
    !      ----------------                                               |
    !      |    |  mat    |                                               |
    !      |    |         | nv1                                           |
    !      ----------------                                               |
    !      |    |         |                                               |
    !      |mat1|         |                                               |
    !      |    |         |                                               |
    !      |    |         | nv2                                           |
    !      ----------------                                               |
    !                                                                     |
    ! For eigenvector-parallelization this needs some communication       |
    ! between the nodes, since this part is created 'column-wise'         |
    ! but needed row-wise.                                                |
    !                                                                     |
    !     n_send(i): number of elements to send to pe #i                  |
    !     n_recv(i): number of elements to receive from pe #i             |
    !     ns_tot,nr_tot : total number of elements to send/receive        |
    !     cs_el,cr_el: send and receive elements                          |
    !     in_pos(xy,n,i): where to put in the n'th element sent by pe #i  |
    !                                                                     |
    ! Based on the mingeselle.F90                                         |
    !                                        Sept. 2020  U.Alekseeva      | 
    !---------------------------------------------------------------------+
#include"./cpp_double.h"
    CLASS(t_mat), INTENT(INOUT) ::mat1
    CLASS(t_mat), INTENT(INOUT) ::mat
    ! ..
    ! .. Local Scalarsxi
    INTEGER n_rank,n_size
    INTEGER ki,kj,kjj,ns_tot,nr_tot,n_p,n_help,i,ii
    INTEGER ns_max,nr_max,np_s,np_r
    INTEGER inext,ifront,req_s,req_r
    INTEGER SUB_COMM
    ! ..
    ! .. Local Arrays
    INTEGER ierr(3)
    INTEGER, ALLOCATABLE :: c_help_size(:,:)
    INTEGER, ALLOCATABLE :: n_send(:),nsr(:),nsc(:)
    INTEGER, ALLOCATABLE :: n_recv(:),n_r(:)
    INTEGER, ALLOCATABLE :: in_pos(:,:,:)
    COMPLEX, ALLOCATABLE :: cs_el(:,:),cr_el(:),b_b(:),c_help(:,:)
    LOGICAL, ALLOCATABLE :: nsl(:)

    INCLUDE 'mpif.h'
    INTEGER stt(MPI_STATUS_SIZE)

    SELECT TYPE (mat1)
    TYPE IS (t_mpimat)
       SELECT TYPE (mat)
       TYPE IS (t_mpimat)
          IF ( (mat1%global_size1 .NE. mat%global_size2) .OR. (mat1%global_size2 .NE. mat%global_size1) ) THEN
             CALL juDFT_error("The matrices do not match",calledby ="mingeselle")
          ENDIF 
          CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, n_rank, ierr)
          CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, n_size, ierr)
          SUB_COMM = mat%blacsdata%mpi_com
          !Set lower part of matrix to zero...
          ii = 0
          DO i = n_rank + 1, MIN(mat%global_size1, mat%global_size2), n_size
            ii = ii + 1
            IF (mat%l_real) THEN
               mat%data_r(i + 1:, ii) = 0.0
               mat1%data_r(i + 1:, ii) = 0.0
               mat1%data_r(i, ii) = 0.0
            ELSE
               mat%data_c(i + 1:, ii) = 0.0
               mat1%data_c(i + 1:, ii) = 0.0
               mat1%data_c(i, ii) = 0.0
            ENDIF
          ENDDO
       END SELECT
    END SELECT
    !
    ! initialize
    !
    ALLOCATE(n_send(0:n_size-1),n_recv(0:n_size-1),n_r(0:n_size-1),nsr(0:n_size-1),nsc(0:n_size-1))
    ALLOCATE(c_help_size(2,0:n_size-1),nsl(0:n_size-1))
    ns_tot = 0
    nr_tot = 0
    n_send = 0
    n_r = 0
    nsc = 0
    c_help_size = 0
    !
    ! determine number of elements to send to other pe's
    ! and calculate the dimensions of c_helpi
    ! rows of c_help correspond to columns of mat1 and vice versa
    !
    DO ki = 1, mat1%matsize2
       kjj = n_rank + 1 + (ki-1)*n_size      ! global column index
       nsr = 0
       nsl = .FALSE.
       DO kj = 1, min(kjj-1,mat1%matsize1)
          ns_tot = ns_tot + 1
          n_p = MOD(kj-1,n_size)
          n_send(n_p) = n_send(n_p) + 1
          nsr(n_p) = nsr(n_p) + 1
          nsl(n_p) = .TRUE.
       ENDDO
       DO n_p = 0,n_size-1
          IF ( c_help_size(2,n_p) < nsr(n_p) ) c_help_size(2,n_p) = nsr(n_p)
          IF ( nsl(n_p)==.TRUE. ) c_help_size(1,n_p) = c_help_size(1,n_p) + 1
       ENDDO
    ENDDO
    !print*, "send", n_rank, ns_tot, n_send
    !
    ! determine number of elements to receive from other pe's
    !
    DO ki = 1, mat%matsize2
       kjj = n_rank + 1 + (ki-1)*n_size      ! global column index
       DO kj = kjj+1, mat%matsize1
          nr_tot = nr_tot + 1
          n_p = MOD(kj-1,n_size)
          n_r(n_p) = n_r(n_p) + 1
       ENDDO
    ENDDO
    !print*, "recv", n_rank, nr_tot, n_r
    !
    ! determine the maximal number of s/r-counts and allocate s/r-arrays
    !
    ns_max = 0
    nr_max = 0
    DO n_p = 0,n_size-1
       ns_max = MAX(ns_max,n_send(n_p))
       nr_max = MAX(nr_max,n_r(n_p))
    ENDDO
    !      WRITE (*,*) ns_max ,nr_max  , n_size, n_rank
    ALLOCATE ( cs_el(ns_max,0:n_size-1),cr_el(nr_max) )
    ALLOCATE ( in_pos(2,nr_max,0:n_size-1) )
    !
    ! sort the elements of the mat1 into the c_help,
    ! resorting them on the way: rows <-> columns
    ! then put them in the send buffers
    !
    ALLOCATE ( c_help(mat1%matsize2,ceiling(real(mat1%matsize1)/n_size)) )
    c_help = cmplx(0.0,0.0)
    DO n_p = 0,n_size-1

       IF (c_help_size(2,n_p) > size(c_help,2))  CALL juDFT_error("allocated c_help is too small",calledby ="mingeselle")
       IF (c_help_size(1,n_p) > size(c_help,1))  CALL juDFT_error("allocated c_help is too small",calledby ="mingeselle")
       !print*, "c_help_size",n_rank, n_p,c_help_size(:,n_p)

       DO ki = 1, c_help_size(1,n_p)
          DO kj = 1, min(ki,c_help_size(2,n_p))
            kjj = (kj-1)*n_size+n_p+1     ! #row of the element in mat1
            IF (n_rank-1 < n_p) THEN 
               c_help(ki,kj) = mat1%data_c(kjj,ki+1)
            ELSE
               c_help(ki,kj) = mat1%data_c(kjj,ki)
            ENDIF
          ENDDO
       ENDDO

       n_help = 0
       DO kj = 1,c_help_size(2,n_p)
          DO ki = kj ,c_help_size(1,n_p)
             n_help = n_help + 1
             cs_el(n_help,n_p) = CONJG(c_help(ki,kj))
          ENDDO
       ENDDO
       IF ( n_help .NE. n_send(n_p)) CALL juDFT_error("Number of elements to send is wrong",calledby ="mingeselle")
   
    ENDDO
    DEALLOCATE ( c_help )
    !CALL juDFT_error("stop",calledby ="mingeselle")
    !
    ! now we look where to put in the received elements
    !
    n_recv = 0
    DO ki = 1, mat%matsize2
       kjj = n_rank + 1 + (ki-1)*n_size      ! global column index
       DO kj = kjj+1, mat%matsize1
          n_p = MOD(kj-1,n_size)
          n_recv(n_p) = n_recv(n_p) + 1
          in_pos(1,n_recv(n_p),n_p) = kj
          in_pos(2,n_recv(n_p),n_p) = ki
       ENDDO
    ENDDO
    DO n_p = 0,n_size-1
       IF (n_recv(n_p)/=n_r(n_p))  CALL juDFT_error("n_recv.NE.n_s" ,calledby ="mingeselle")
    ENDDO
    !
    ! Mandaliet, mandaliet, min geselle kumme niet
    !
    ifront = ibefore(n_size,n_rank)
    inext  = iafter (n_size,n_rank)
    DO n_p = 0,n_size-1
       !
       ! determine pe's to send to and to receive from
       !
       np_s = MOD(inext +n_p,n_size)
       np_r = MOD(ifront-n_p,n_size)
       IF (np_r.LT.0) np_r = np_r + n_size
       !
       ! send section: local rows i with mod(i-1,np) = np_s will be sent to proc np_s
       !

       IF (np_s.NE.n_rank) THEN
          CALL MPI_ISEND(cs_el(1,np_s),n_send(np_s), CPP_MPI_COMPLEX,&
               np_s,n_rank,SUB_COMM,req_s,ierr)
          !          write (*,*) n_rank,'sends',n_send(np_s),'to',np_s
          !          write (*,'(i2,10f10.7)') n_rank,(real(cs_el(ki,np_s)),ki=1,10)
       ENDIF

       !
       ! receive section : local rows i  with mod(i-1,np) = np_r will be received from np_r
       ! ... skipped, if update matrix from local data:
       !
       IF (np_r.NE.n_rank) THEN
          CALL MPI_IRECV(cr_el,n_recv(np_r),CPP_MPI_COMPLEX, MPI_ANY_SOURCE,np_r,SUB_COMM,req_r,ierr)
          CALL MPI_WAIT(req_s,stt,ierr)
          CALL MPI_WAIT(req_r,stt,ierr)
          !          write (*,*) n_rank,'recvs',ierr,n_p,np_r
          !          write(*,*) n_rank,'receives',n_recv(np_r),'from',np_r
          !          write (*,'(i2,10f10.7)') n_rank,(real(cr_el(ki)),ki=1,10)
          !
          ! now update the matrix aa()
          !
          !IF (l_real) THEN
          !   aa_r(in_pos(:n_recv(np_r),np_r)) = aa_r(in_pos(:n_recv(np_r),np_r)) + cr_el(:n_recv(np_r))
          !ELSE
             !aa_c(in_pos(:n_recv(np_r),np_r)) = aa_c(in_pos(:n_recv(np_r),np_r)) + cr_el(:n_recv(np_r))
             DO ki = 1,n_recv(np_r)
                mat%data_c(in_pos(1,ki,np_r),in_pos(2,ki,np_r)) = cr_el(ki)
             ENDDO
          !ENDIF
       ELSE
          !IF (l_real) THEN
          !   aa_r(in_pos(:n_recv(np_r),np_r)) = aa_r(in_pos(:n_recv(np_r),np_r)) + cs_el(:n_recv(np_r),np_s)
          !ELSE
             !aa_c(in_pos(:n_recv(np_r),np_r)) = aa_c(in_pos(:n_recv(np_r),np_r)) + cs_el(:n_recv(np_r),np_s)
             DO ki = 1,n_recv(np_r)
                mat%data_c(in_pos(1,ki,np_r),in_pos(2,ki,np_r)) = cs_el(ki,np_s)
             ENDDO
          !ENDIF
       ENDIF
       !         CALL MPI_BARRIER(SUB_COMM,ierr)
    ENDDO

    DEALLOCATE (cs_el,cr_el,in_pos)

  END SUBROUTINE mingeselle
  !
  !-------------------------------------------------------------
  !
  INTEGER FUNCTION ibefore(np, p)
    !
    ! Determine (in a ring structure) which is the front process
    !
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: np  !  number of processes
    INTEGER, INTENT (IN) :: p   !  current processes

    IF ( p > 0 ) THEN
       ibefore = p-1
    ELSE
       ibefore = np-1
    ENDIF

  END FUNCTION ibefore
  !
  !-------------------------------------------------------------
  !
  INTEGER FUNCTION iafter(np, p)
    !
    ! Determine (in a ring structure) which is the next process
    !
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: np  !  number of processes
    INTEGER, INTENT (IN) :: p   !  current processes

    IF ( p < np-1 ) THEN
       iafter = p+1
    ELSE
       iafter = 0
    ENDIF

  END FUNCTION iafter
  !
  !-------------------------------------------------------------
  !
END MODULE m_mingeselle
