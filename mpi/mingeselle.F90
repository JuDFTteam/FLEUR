MODULE m_mingeselle
  USE m_juDFT
CONTAINS
  SUBROUTINE mingeselle(SUB_COMM,n_size,n_rank,nv, ahelp,l_real,aa_r,aa_c)
    !------------------------------------------------------------------+
    !                                                                  |
    ! Transfers the spin-down/spin-up part , upper triangle of the     |
    ! MT-hamiltonian from the help-array ``ahelp'' to the H-matrix.    |
    ! For eigenvector-parallelization this needs some communication    |
    ! between the nodes, since this part is created 'column-wise'      |
    ! but needed row-wise.                                             |
    !                                                                  |
    !     n_s(i): number of elements to send to pe #i                  |
    !     n_r(i): number of elements to receive from pe #i             |
    !     ns_tot,nr_tot : total number of elements to send/receive     |
    !     n2_start: pe that has first column of 2nd spin part          |
    !     cs_el,cr_el: send and receive elements                       |
    !     in_pos(n,i): where to put in the n'th element sent by pe #i  |
    !                                                                  |
    !------------------------------------------------------------------+
#include"./cpp_double.h"

    IMPLICIT NONE
    ! ..
    ! .. Scalar Arguments
    INTEGER, INTENT (IN) :: n_size,n_rank,SUB_COMM
    ! ..
    ! .. Array Arguments
    INTEGER, INTENT (IN)    :: nv(2)
    COMPLEX, INTENT (INOUT) :: ahelp(:)!(m_ahelp)
    LOGICAL, INTENT (IN)    :: l_real
    REAL,   OPTIONAL, INTENT (INOUT) :: aa_r(:)!(matsize)
    COMPLEX,OPTIONAL, INTENT (INOUT) :: aa_c(:)
    ! ..
    ! .. Local Scalars
    INTEGER ki,kj,ns_tot,nr_tot,n_p,n2_start,n_help
    INTEGER ns_max,nr_max,n_pos,np_s,np_r,nv_s,ii,i
    INTEGER inext,ifront,req_s,req_r
    ! ..
    ! .. Local Arrays
    INTEGER n_s(0:n_size-1),n_send(0:n_size-1)
    INTEGER n_r(0:n_size-1),n_recv(0:n_size-1),ierr(3)
    INTEGER, ALLOCATABLE :: in_pos(:,:)
    COMPLEX, ALLOCATABLE :: cs_el(:,:),cr_el(:),b_b(:),c_help(:,:)

    INCLUDE 'mpif.h'
    INTEGER stt(MPI_STATUS_SIZE)
    ! ..
    !
    ! kick out the diagonal elements of ahelp
    !
    i  = 0
    ii = 0
    DO ki =  n_rank+1, nv(1), n_size
       DO kj = 1,ki - 1
          i  =  i + 1
          ii = ii + 1
          ahelp(i) = ahelp(ii)
       END DO
       ii = ii + 1
    ENDDO
    !
    ! initialize
    !
    ns_tot = 0
    nr_tot = 0
    DO n_p = 0,n_size-1
       n_s(n_p) = 0
       n_r(n_p) = 0
    ENDDO
    !
    ! determine number of elements to send to other pe's
    !
    n2_start = MOD(nv(1),n_size) - 1
    DO ki = 1, nv(1)
       IF ( MOD(ki-1,n_size).EQ.n_rank ) THEN
          DO kj = 1, ki-1
             ns_tot = ns_tot + 1
             n_p = MOD((kj+n2_start),n_size)
             n_s(n_p) = n_s(n_p) + 1
          ENDDO
       ENDIF
    ENDDO
    !
    ! determine number of elements to receive from other pe's
    !
    DO ki = 1, nv(2)
       IF ( MOD(ki+nv(1)-1,n_size).EQ.n_rank ) THEN
          DO kj = ki+1, nv(2)
             nr_tot = nr_tot + 1
             n_p = MOD(kj-1,n_size)
             n_r(n_p) = n_r(n_p) + 1
          ENDDO
       ENDIF
    ENDDO
    !
    !      WRITE (*,*) ns_tot,(n_s(n_p),n_p=0,n_size-1)
    !      WRITE (*,*) nr_tot,(n_r(n_p),n_p=0,n_size-1)
    !
    ! determine the maximal number of s/r-counts and allocate s/r-arrays
    !
    ns_max = 0
    nr_max = 0
    DO n_p = 0,n_size-1
       ns_max = MAX(ns_max,n_s(n_p))
       nr_max = MAX(nr_max,n_r(n_p))
    ENDDO
    !      WRITE (*,*) ns_max ,nr_max  , n_size, n_rank
    ALLOCATE ( cs_el(ns_max,0:n_size-1),cr_el(nr_max), in_pos(nr_max,0:n_size-1) )
    !
    ! sort the elements of aahelp-array into the send-arrays
    !
    n_help = 0
    DO n_p = 0,n_size-1
       n_send(n_p) = 0
    ENDDO
    DO ki = 1, nv(1)
       IF ( MOD(ki-1,n_size).EQ.n_rank ) THEN
          DO kj = 1, ki-1
             n_help = n_help + 1
             n_p = MOD((kj+n2_start),n_size)
             n_send(n_p) = n_send(n_p) + 1
             cs_el(n_send(n_p),n_p) = ahelp(n_help)
          ENDDO
       ENDIF
    ENDDO
    IF (n_help/=ns_tot)  CALL juDFT_error("n_help.NE.ns_to         t",calledby ="mingeselle")
    DO n_p = 0,n_size-1
       IF (n_send(n_p)/=n_s(n_p))  CALL juDFT_error("n_send.NE.n_s" ,calledby ="mingeselle")
    ENDDO
    !
    ! resort send array: rows <-> columns
    !
    DO n_p = 0,n_size-1
       nv_s = NINT(SQRT(2.0*n_send(n_p))-0.5)
       ALLOCATE ( c_help(nv_s,nv_s) )

       n_help = 0
       DO ki = 1,nv_s
          DO kj = 1,ki
             n_help = n_help + 1
             c_help(ki,kj) = cs_el(n_help,n_p)
          ENDDO
       ENDDO

       n_help = 0
       DO kj = 1,nv_s
          DO ki = kj ,nv_s
             n_help = n_help + 1
             cs_el(n_help,n_p) = c_help(ki,kj)
          ENDDO
       ENDDO

       DEALLOCATE ( c_help )
    ENDDO
    !
    ! now we look where to put in the received elements
    !
    n_pos = 0
    DO n_p = 0,n_size-1
       n_recv(n_p) = 0
    ENDDO
    DO ki = 1, nv(1)+nv(2)
       IF ( MOD(ki-1,n_size).EQ.n_rank ) THEN
          DO kj = 1, ki
             n_pos = n_pos + 1 
             IF ( ki.GT.nv(1) ) THEN
                IF ((kj.GT.ki-nv(1)).AND.(kj.LE.nv(1))) THEN
                   n_p = MOD(kj-1,n_size)
                   n_recv(n_p) = n_recv(n_p) + 1
                   in_pos(n_recv(n_p),n_p) = n_pos
                ENDIF
             ENDIF
          ENDDO
       ENDIF
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
          IF (l_real) THEN
             aa_r(in_pos(:n_recv(np_r),np_r)) = aa_r(in_pos(:n_recv(np_r),np_r)) + cr_el(:n_recv(np_r))
          ELSE
             aa_c(in_pos(:n_recv(np_r),np_r)) = aa_c(in_pos(:n_recv(np_r),np_r)) + cr_el(:n_recv(np_r))
          ENDIF
       ELSE
          IF (l_real) THEN
             aa_r(in_pos(:n_recv(np_r),np_r)) = aa_r(in_pos(:n_recv(np_r),np_r)) + cs_el(:n_recv(np_r),np_s)
          ELSE
             aa_c(in_pos(:n_recv(np_r),np_r)) = aa_c(in_pos(:n_recv(np_r),np_r)) + cs_el(:n_recv(np_r),np_s)
          ENDIF
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
