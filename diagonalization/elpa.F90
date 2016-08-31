!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

#ifdef CPP_INVERSION
#define CPP_CHOLESKY cholesky_real
#define CPP_invert_trm invert_trm_real
#define CPP_solve_evp solve_evp_real
#define CPP_solve_evp_2stage solve_evp_real_2stage
#define CPP_REALCOMPLEX real
#define CPP_transpose pdtran
#define CPP_ONE 1.0
#define CPP_ZERO 0.0
#define CPP_mult mult_at_b_real
#else
#define CPP_CHOLESKY cholesky_complex
#define CPP_invert_trm invert_trm_complex
#define CPP_solve_evp solve_evp_complex
#define CPP_solve_evp_2stage solve_evp_complex_2stage
#define CPP_REALCOMPLEX complex
#define CPP_transpose pztranc
#define CPP_ONE cmplx(1.,0.)
#define CPP_ZERO cmplx(0.,0.)
#define CPP_mult mult_ah_b_complex
#endif
MODULE m_elpa
  PRIVATE
  !Module to call elpa library for parallel diagonalization
  !uses ssubredist1/2 for redistribution

  PUBLIC elpa

CONTAINS

  SUBROUTINE elpa(m,n, SUB_COMM, a,b, z,eig,num)
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
    USE m_juDFT
    USE elpa1
#ifdef CPP_ELPA2
    USE elpa2
#endif
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER, INTENT (IN)                  :: m,n
    INTEGER, INTENT (IN)                  :: SUB_COMM
    INTEGER, INTENT (INOUT)               :: num
    REAL,    INTENT   (OUT)               :: eig(:)
    CPP_REALCOMPLEX, ALLOCATABLE, INTENT (INOUT)  :: a(:),b(:)
    CPP_REALCOMPLEX, INTENT   (OUT)               :: z(:,:)
    !
    !...  Local variables
    !
    INTEGER           :: sc_desc(9) !SCALAPACK DESCRIPTOR
    INTEGER           :: np,nb,myid,npcol,nprow,mycol,myrow,myrowssca
    INTEGER           :: mycolssca,ictextblacs,mpi_comm_rows
    INTEGER           :: mpi_comm_cols
    INTEGER           :: num2

    LOGICAL           :: ok
    INTEGER           :: err,ierr
    INTEGER           :: i,k,n_col,n_row
    ! large distributed Matrices
    REAL,ALLOCATABLE     :: eig2(:)
    CPP_REALCOMPLEX, ALLOCATABLE :: asca(:,:), bsca(:,:),tmp2(:,:)
    CPP_REALCOMPLEX, ALLOCATABLE :: eigvec(:,:)
    INTEGER, EXTERNAL :: numroc, indxl2g  !SCALAPACK functions


    CALL priv_create_blacsgrid(sub_comm,m, np,nb,myid,npcol,nprow,mycol,myrow,myrowssca,mycolssca,&
         ictextblacs,sc_desc,mpi_comm_rows,mpi_comm_cols)

    !     Number of columns the local process gets in ScaLAPACK distribution
    ALLOCATE ( asca(myrowssca,mycolssca), stat=err )
    IF (err/=0) CALL juDFT_error("allocate asca", calledby="chani_elpa")
    asca=0.0

    ALLOCATE ( bsca(myrowssca,mycolssca), stat=err  )
    IF (err/=0) CALL juDFT_error("allocate bsca", calledby="chani_elpa")
    bsca=0.0
    !
    ! transfer to ScaLAPACK block cyclic 2d distribution with
    ! block-size nb
    !
    !print *,"Before subredist"
    CALL subredist1(a,m,asca,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb )
    CALL subredist1(b,m,bsca,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb )

    ! for saving memory one might deallocate(a,b)
    !print *,"Before Allocate"
    ALLOCATE ( eig2(m), stat=err ) ! The eigenvalue array for ScaLAPACK
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "eig2"', calledby ='elpa')


    ALLOCATE ( eigvec(myrowssca,mycolssca), stat=err ) ! Eigenvectors for ScaLAPACK
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "eigvec"',calledby ='elpa')

    ALLOCATE ( tmp2(myrowssca,mycolssca), stat=err ) ! tmp_array
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "tmp2"', calledby ='elpa')
    !print *,"Before elpa"

    !ELPA -start here
    ! Solve generalized problem
    !
    ! 1. Calculate Cholesky factorization of Matrix B = U**T * U
    !    and invert triangular matrix U
    !
    ! Please note: cholesky_complex/invert_trm_complex are not trimmed for speed.
    ! The only reason having them is that the Scalapack counterpart
    ! PDPOTRF very often fails on higher processor numbers for unknown reasons!
#ifdef CPP_ELPA_NEW
    CALL CPP_CHOLESKY (m,bsca,SIZE(bsca,1),nb,mycolssca,mpi_comm_rows,mpi_comm_cols,.false.,ok)
    CALL CPP_invert_trm(m,bsca,SIZE(bsca,1),nb,mycolssca,mpi_comm_rows,mpi_comm_cols,.false.,ok)
#else
    CALL CPP_CHOLESKY (m,bsca,SIZE(bsca,1),nb, mpi_comm_rows,mpi_comm_cols)
    CALL CPP_invert_trm(m, bsca, SIZE(bsca,1), nb, mpi_comm_rows, mpi_comm_cols)
#endif

    ! 2. Calculate U**-T * A * U**-1

    ! 2a. eigvec = U**-T * A

    ! A is only set in the upper half, solve_evp_real needs a full matrix
    ! Set lower half from upper half
    CALL CPP_transpose (m,m,CPP_ONE,asca,1,1,sc_desc, CPP_ZERO,eigvec,1,1,sc_desc)
    DO i=1,SIZE(asca,1)
       ! Get global column corresponding to i and number of local rows up to
       ! and including the diagonal, these are unchanged in A
       n_col = indxl2g(i,     nb, mycol, 0, npcol)
       n_row = numroc (n_col, nb, myrow, 0, nprow)
       asca(n_row+1:myrowssca,i) = eigvec(n_row+1:myrowssca,i)
    ENDDO

    CALL CPP_mult ('U', 'L',m, m,bsca,myrowssca,asca,SIZE(asca,1),nb, mpi_comm_rows, mpi_comm_cols,eigvec,myrowssca)

    ! 2b. tmp2 = eigvec**T
    CALL CPP_transpose(m,m,CPP_ONE,eigvec,1,1,sc_desc,CPP_ZERO,tmp2,1,1,sc_desc)

    ! 2c. A =  U**-T * tmp2 ( = U**-T * Aorig * U**-1 )
    CALL CPP_mult ('U', 'U', m, m, bsca, SIZE(bsca,1), tmp2,&
         SIZE(tmp2,1),nb, mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1))

    ! A is only set in the upper half, solve_evp_real needs a full matrix
    ! Set lower half from upper half

    CALL CPP_transpose(m,m,CPP_ONE,asca,1,1,sc_desc, CPP_ZERO,eigvec,1,1,sc_desc)


    DO i=1,SIZE(asca,1)
       ! Get global column corresponding to i and number of local rows up to
       ! and including the diagonal, these are unchanged in A
       n_col = indxl2g(i,     nb, mycol, 0, npcol)
       n_row = numroc (n_col, nb, myrow, 0, nprow)
       asca(n_row+1:myrowssca,i) = eigvec(n_row+1:myrowssca,i)
    ENDDO

    ! 3. Calculate eigenvalues/eigenvectors of U**-T * A * U**-1
    !    Eigenvectors go to eigvec
    num2=num
#ifdef CPP_ELPA_NEW
#ifdef CPP_ELPA2
    err=CPP_solve_evp_2stage(m,num2,asca,SIZE(asca,1),&
         eig2,eigvec,SIZE(asca,1), nb,mycolssca, mpi_comm_rows, mpi_comm_cols,sub_comm)
#else
    err=CPP_solve_evp(m, num2,asca,SIZE(asca,1),&
         eig2,eigvec, SIZE(asca,1), nb,mpi_comm_rows, mpi_comm_cols)
#endif
#else
#ifdef CPP_ELPA2
    CALL CPP_solve_evp_2stage(m,&
         num2,asca,SIZE(asca,1),eig2,eigvec,SIZE(asca,1), nb,mpi_comm_rows, mpi_comm_cols,sub_comm)
#else
    CALL CPP_solve_evp(m, num2,&
         asca,SIZE(asca,1),eig2,eigvec, SIZE(asca,1), nb,mpi_comm_rows, mpi_comm_cols)
#endif
#endif

    ! 4. Backtransform eigenvectors: Z = U**-1 * eigvec

    ! mult_ah_b_complex needs the transpose of U**-1, thus tmp2 = (U**-1)**T
    CALL CPP_transpose(m,m,CPP_ONE,bsca,1,1,sc_desc,CPP_ZERO,tmp2,1,1,sc_desc)

    CALL CPP_mult ('L', 'N',m, num2, tmp2, SIZE(asca,1),&
         eigvec, SIZE(asca,1),nb,mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1))

    ! END of ELPA stuff
    CALL BLACS_GRIDEXIT(ictextblacs,ierr)
    CALL MPI_COMM_FREE(mpi_comm_rows,ierr)
    CALL MPI_COMM_FREE(mpi_comm_cols,ierr)
    !print *,"elpa done"

    !
    !     Put those eigenvalues expected by chani to eig, i.e. for
    !     process i these are eigenvalues i+1, np+i+1, 2*np+i+1...
    !     Only num=num2/np eigenvalues per process
    !
    num=FLOOR(REAL(num2)/np)
    IF (myid.LT.num2-(num2/np)*np) num=num+1
    k=1
    DO i=myid+1,num2,np
       eig(k)=eig2(i)
       k=k+1
    ENDDO

    !write(*,"(a,i5,99f10.4)") "ELPA:",myid,eig

    !
    !     Redistribute eigvec from ScaLAPACK distribution to each process
    !     having all eigenvectors corresponding to his eigenvalues as above
    !
    CALL subredist2(z,m,num2,asca,myrowssca,SUB_COMM,nprow,npcol, myid,ierr,nb)
    !

    DEALLOCATE ( asca )
    DEALLOCATE ( bsca )
    DEALLOCATE (eigvec, tmp2, eig2)

  END SUBROUTINE elpa

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
#ifdef CPP_ELPA_NEW
    ierr=get_elpa_row_col_comms(mpi_subcom, myrowblacs, mycolblacs,mpi_comm_rows, mpi_comm_cols)
#else
    CALL get_elpa_row_col_comms(mpi_subcom, myrowblacs, mycolblacs,mpi_comm_rows, mpi_comm_cols)
#endif
    !print *,"creating ELPA comms  --  done"

    !Create the descriptors
    CALL descinit(sc_desc,m,m,nb,nb,0,0,ictextblacs,myrowssca,ierr)
    IF (ierr /=0 ) CALL juDFT_error('descinit1 failed',calledby='elpa')

  END SUBROUTINE priv_create_blacsgrid

END MODULE m_elpa
