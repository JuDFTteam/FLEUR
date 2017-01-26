!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_elpa
  PRIVATE
#ifdef CPP_ELPA  
  interface elpa_diag
     module procedure elpa_r,elpa_c
  end interface elpa_diag
  
  !Module to call elpa library for parallel diagonalization
  !uses ssubredist1/2 for redistribution

  PUBLIC elpa_diag

CONTAINS
  ! First the real version of the code
#define CPP_CHOLESKY cholesky_real
#define CPP_invert_trm invert_trm_real
#define CPP_solve_evp_1stage solve_evp_real_1stage
#define CPP_solve_evp solve_evp_real
#define CPP_solve_evp_2stage solve_evp_real_2stage
#define CPP_REALCOMPLEX real
#define CPP_transpose pdtran
#define CPP_ONE 1.0
#define CPP_ZERO 0.0
#define CPP_mult mult_at_b_real
#define CPP_REAL
  SUBROUTINE elpa_r(m, SUB_COMM, a,b, z,eig,num)
    ! 
    !----------------------------------------------------
    !- Parallel eigensystem solver - driver routine based on chani; dw'12
    !
    ! m ........ actual (=leading) dimension of full a & b matrices
    !            must be problem size, as input a, b  are one-dimensional
    !            and shall be redistributed to two-dimensional matrices
    !            actual (=leading) dimension of eigenvector z(,)
    ! n ........ number of columns of full (sub)matrix ( about n/np)
    ! SUB_COMM.. communicator for MPI
    ! a,b   .... packed (sub)matrices, here expanded to non-packed
    ! z,eig .... eigenvectors and values, output
    ! num ...... number of ev's searched (and found) on this node
    !            On input, overall number of ev's searched,
    !            On output, local number of ev's found
    !
    !----------------------------------------------------
    !
  
#include "elpa_cpp.F90"

  END SUBROUTINE elpa_r


  ! Now the complex version
#define CPP_CHOLESKY cholesky_complex
#define CPP_invert_trm invert_trm_complex
#define CPP_solve_evp_1stage solve_evp_complex_1stage
#define CPP_solve_evp solve_evp_complex
#define CPP_solve_evp_2stage solve_evp_complex_2stage
#define CPP_REALCOMPLEX complex
#define CPP_transpose pztranc
#define CPP_ONE cmplx(1.,0.)
#define CPP_ZERO cmplx(0.,0.)
#define CPP_mult mult_ah_b_complex
#undef CPP_REAL
 SUBROUTINE elpa_c(m, SUB_COMM, a,b, z,eig,num)
    ! 
    !----------------------------------------------------
    !- Parallel eigensystem solver - driver routine based on chani; dw'12
    !
    ! m ........ actual (=leading) dimension of full a & b matrices
    !            must be problem size, as input a, b  are one-dimensional
    !            and shall be redistributed to two-dimensional matrices
    !            actual (=leading) dimension of eigenvector z(,)
    ! n ........ number of columns of full (sub)matrix ( about n/np)
    ! SUB_COMM.. communicator for MPI
    ! a,b   .... packed (sub)matrices, here expanded to non-packed
    ! z,eig .... eigenvectors and values, output
    ! num ...... number of ev's searched (and found) on this node
    !            On input, overall number of ev's searched,
    !            On output, local number of ev's found
    !
    !----------------------------------------------------
    !
  
#include "elpa_cpp.F90"

  END SUBROUTINE elpa_c

  SUBROUTINE priv_create_blacsgrid(mpi_subcom,m, np,nb,myid,npcol,nprow,mycol,&
       myrow,myrowssca,mycolssca, ictextblacs,sc_desc,mpi_comm_rows,mpi_comm_cols)
    USE m_juDFT
    USE elpa1
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,INTENT(IN) :: mpi_subcom
    INTEGER,INTENT(IN) :: m
    INTEGER,INTENT(OUT):: np,nb, npcol,nprow,myrowssca,mycolssca
    INTEGER,INTENT(OUT):: myrow,mycol,myid
    INTEGER,INTENT(OUT):: ictextblacs,sc_desc(:)
    INTEGER,INTENT(OUT):: mpi_comm_rows,mpi_comm_cols


    INTEGER     :: iamblacs,npblacs
    INTEGER     :: nprow2,npcol2,myrowblacs,mycolblacs
    INTEGER     :: k,i,j
    INTEGER     :: ierr

    INTEGER,ALLOCATABLE :: iblacsnums(:),ihelp(:),iusermap(:,:)

    EXTERNAL descinit, blacs_get
    EXTERNAL blacs_pinfo, blacs_gridinit

    !Determine rank and no of processors
    CALL MPI_COMM_RANK(mpi_subcom,myid,ierr)
    CALL MPI_COMM_SIZE(mpi_subcom,np,ierr)

    !print *,"priv_create_blacsgrid"
    ! determine block size
    !
    nb = 20
    IF (m.GT.2048)   nb = 30 !2
    IF (m.GT.8*2048) nb = 60 !4

    ! compute processor grid, as square as possible
    ! If not square with more rows than columns

    distloop: DO j=INT(SQRT(REAL(np))),1,-1
       IF ( (np/j) * j == np) THEN
          npcol = np/j
          nprow = j
          EXIT distloop
       ENDIF
    ENDDO distloop
    ALLOCATE(iblacsnums(np),ihelp(np),iusermap(nprow,npcol))

    !   An nprow*npcol processor grid will be created
    !   Row and column index myrow, mycol of this processor in the grid
    !   and distribution of A and B in ScaLAPACK
    !   The local processor will get myrowssca rows and mycolssca columns
    !   of A and B
    !

    myrow = myid/npcol  ! my row number in the BLACS nprow*npcol grid
    mycol = myid -(myid/npcol)*npcol  ! my column number in the BLACS nprow*npcol grid
    !
    !  Now allocate Asca to put the elements of Achi or receivebuffer to
    !
    myrowssca=(m-1)/(nb*nprow)*nb+ MIN(MAX(m-(m-1)/(nb*nprow)*nb*nprow-nb*myrow,0),nb)
    !     Number of rows the local process gets in ScaLAPACK distribution
    mycolssca=(m-1)/(nb*npcol)*nb+ MIN(MAX(m-(m-1)/(nb*npcol)*nb*npcol-nb*mycol,0),nb)



    !Get BLACS ranks for all MPI ranks
    CALL BLACS_PINFO(iamblacs,npblacs)  ! iamblacs = local process rank (e.g. myid)
    ! npblacs  = number of available processes
    iblacsnums=-2
    ihelp=-2
    ihelp(myid+1)=iamblacs ! Get the Blacs id corresponding to the MPI id
    !print *,"ALLREDUCE:",mpi_subcom
    CALL MPI_ALLREDUCE(ihelp, iblacsnums, np,MPI_INTEGER,MPI_MAX,mpi_subcom,ierr)
    IF (ierr.NE.0)CALL juDFT_error('Error in allreduce for BLACS nums' ,calledby='elpa')

    !     iblacsnums(i) is the BLACS-process number of MPI-process i-1
    k = 1
    DO i = 1, nprow
       DO j = 1, npcol
          iusermap(i,j) = iblacsnums(k)
          k = k + 1
       ENDDO
    ENDDO
    !Get the Blacs default context
    CALL BLACS_GET(0,0,ictextblacs)
    ! Create the Grid
    CALL BLACS_GRIDMAP(ictextblacs,iusermap,nprow,nprow,npcol)
    !     Now control, whether the BLACS grid is the one we wanted
    CALL BLACS_GRIDINFO(ictextblacs, nprow2,npcol2,myrowblacs,mycolblacs)
    IF (nprow2 /= nprow) THEN
       WRITE(6,*) 'Wrong number of rows in BLACS grid'
       WRITE(6,*) 'nprow=',nprow,' nprow2=',nprow2
       CALL juDFT_error('Wrong number of rows in BLACS grid',calledby= 'elpa')
    ENDIF
    IF (npcol2 /= npcol) THEN
       WRITE(6,*) 'Wrong number of columns in BLACS grid'
       WRITE(6,*) 'npcol=',npcol,' npcol2=',npcol2
       CALL juDFT_error('Wrong number of columns in BLACS grid', calledby ='elpa')

    ENDIF

    !Create communicators for ELPA
    !print *,"creating ELPA comms"
#if (defined CPP_ELPA_201605003)||defined(CPP_ELPA_NEW)
    ierr=get_elpa_row_col_comms(mpi_subcom, myrowblacs, mycolblacs,mpi_comm_rows, mpi_comm_cols)
#else
    CALL get_elpa_row_col_comms(mpi_subcom, myrowblacs, mycolblacs,mpi_comm_rows, mpi_comm_cols)
#endif
    !print *,"creating ELPA comms  --  done"

    !Create the descriptors
    CALL descinit(sc_desc,m,m,nb,nb,0,0,ictextblacs,myrowssca,ierr)
    IF (ierr /=0 ) CALL juDFT_error('descinit1 failed',calledby='elpa')

  END SUBROUTINE priv_create_blacsgrid
#else
  LOGICAL :: elpa_used=.false.
#endif
end MODULE m_elpa

