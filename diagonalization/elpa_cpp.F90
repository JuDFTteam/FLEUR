!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
    USE m_juDFT
    USE m_subredist1
    USE m_subredist2
    USE elpa1
#ifdef CPP_ELPA2
    USE elpa2
#endif
#ifdef CPP_ELPA_201705003    
    USE elpa
#endif
IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER, INTENT (IN)                  :: m
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
#ifdef CPP_ELPA_201705003    
    INTEGER :: kernel
    CLASS(elpa_t),pointer :: elpa_obj

    err = elpa_init(20170403)
    elpa_obj => elpa_allocate()
#endif

    num2=num
    CALL priv_create_blacsgrid(sub_comm,m, np,nb,myid,npcol,nprow,mycol,myrow,myrowssca,mycolssca,&
         ictextblacs,sc_desc,mpi_comm_rows,mpi_comm_cols)
#ifdef CPP_ELPA_201705003    
    CALL elpa_obj%set("na", m, err)
    CALL elpa_obj%set("nev", num2, err)
    CALL elpa_obj%set("local_nrows", myrowssca, err)
    CALL elpa_obj%set("local_ncols", mycolssca, err)
    CALL elpa_obj%set("nblk", nb, err)
    CALL elpa_obj%set("mpi_comm_parent", sub_comm, err)
    CALL elpa_obj%set("process_row", myrow, err)
    CALL elpa_obj%set("process_col", mycol, err)
#ifdef CPP_ELPA2
    CALL elpa_obj%set("solver", ELPA_SOLVER_2STAGE)
#else
    CALL elpa_obj%set("solver", ELPA_SOLVER_1STAGE)
#endif
    err = elpa_obj%setup()
    if (myid == 0) then
       call elpa_obj%get("complex_kernel", kernel)
       print *, "elpa uses " // elpa_int_value_to_string("complex_kernel", kernel) // " kernel"
    endif
#endif


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
#ifdef CPP_REAL
    CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb, achi_r=a, asca_r=asca )
    CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb, achi_r=b, asca_r=bsca)
#else
    CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb, achi_c=a, asca_c=asca )
    CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb, achi_c=b, asca_c=bsca)
#endif
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

#if defined(CPP_ELPA_201705003)
    CALL elpa_obj%cholesky(bsca, err)
    CALL elpa_obj%invert_triangular(bsca, err)
#elif defined(CPP_ELPA_201605003) || defined(CPP_ELPA_201605004)
    ok=CPP_CHOLESKY (m,bsca,SIZE(bsca,1),nb,mycolssca,mpi_comm_rows,mpi_comm_cols,.false.)
    ok=CPP_invert_trm(m,bsca,SIZE(bsca,1),nb,mycolssca,mpi_comm_rows,mpi_comm_cols,.false.)
#elif defined CPP_ELPA_NEW
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

#if defined (CPP_ELPA_201705003)
    CALL elpa_obj%hermitian_multiply('U','L', m, bsca, asca, myrowssca, mycolssca, eigvec, myrowssca, mycolssca, err)
#elif defined (CPP_ELPA_201605004)
    ok=CPP_mult ('U', 'L',m, m,bsca,myrowssca,mycolssca,asca,SIZE(asca,1),SIZE(asca,2),nb,&
         mpi_comm_rows, mpi_comm_cols,eigvec,myrowssca,mycolssca)
#elif defined (CPP_ELPA_201605003)
    ok=CPP_mult ('U', 'L',m, m,bsca,myrowssca,asca,SIZE(asca,1),nb, mpi_comm_rows, mpi_comm_cols,eigvec,myrowssca)
#else
    CALL CPP_mult ('U', 'L',m, m,bsca,myrowssca,asca,SIZE(asca,1),nb, mpi_comm_rows, mpi_comm_cols,eigvec,myrowssca)
#endif

    ! 2b. tmp2 = eigvec**T
    CALL CPP_transpose(m,m,CPP_ONE,eigvec,1,1,sc_desc,CPP_ZERO,tmp2,1,1,sc_desc)

    ! 2c. A =  U**-T * tmp2 ( = U**-T * Aorig * U**-1 )
#if defined (CPP_ELPA_201705003)
    CALL elpa_obj%hermitian_multiply('U','U', m, bsca, tmp2, myrowssca, mycolssca, asca, myrowssca, mycolssca, err)
#elif defined (CPP_ELPA_201605004)
    ok=CPP_mult ('U', 'U', m, m, bsca, SIZE(bsca,1),SIZE(bsca,2), tmp2,&
         SIZE(tmp2,1),SIZE(tmp2,2),nb, mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1),SIZE(asca,2))
#elif defined (CPP_ELPA_201605003)
    ok=CPP_mult ('U', 'U', m, m, bsca, SIZE(bsca,1), tmp2,&
         SIZE(tmp2,1),nb, mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1))
#else
    CALL CPP_mult ('U', 'U', m, m, bsca, SIZE(bsca,1), tmp2,&
         SIZE(tmp2,1),nb, mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1))
#endif

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
#if defined (CPP_ELPA_201705003)
    CALL elpa_obj%eigenvectors(asca, eig2, eigvec, err)
#elif defined(CPP_ELPA_201605003) || defined(CPP_ELPA_201605004)
#ifdef CPP_ELPA2
    ok=CPP_solve_evp_2stage(m,num2,asca,SIZE(asca,1),&
         eig2,eigvec,SIZE(asca,1), nb,mycolssca, mpi_comm_rows, mpi_comm_cols,sub_comm)
#else
    ok=CPP_solve_evp_1stage(m, num2,asca,SIZE(asca,1),&
         eig2,eigvec, SIZE(asca,1), nb,mycolssca,mpi_comm_rows, mpi_comm_cols)
#endif
#elif defined CPP_ELPA_NEW
#ifdef CPP_ELPA2
    err=CPP_solve_evp_2stage(m,num2,asca,SIZE(asca,1),&
         eig2,eigvec, SIZE(asca,1), nb,mycolssca, mpi_comm_rows, mpi_comm_cols,sub_comm)
#else
    err=CPP_solve_evp(m, num2,asca,SIZE(asca,1),&
         eig2,eigvec, SIZE(asca,1), nb,mycolssca, mpi_comm_rows, mpi_comm_cols)
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

#if defined (CPP_ELPA_201705003)
    CALL elpa_obj%hermitian_multiply('L','N', num2, tmp2, eigvec, myrowssca, mycolssca, asca, myrowssca, mycolssca, err)
#elif defined (CPP_ELPA_201605004)
    ok= CPP_mult ('L', 'N',m, num2, tmp2, SIZE(asca,1),SIZE(asca,2),&
         eigvec, SIZE(asca,1),SIZE(asca,2),nb,mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1),SIZE(asca,2))
#elif CPP_ELPA_201605003
    ok= CPP_mult ('L', 'N',m, num2, tmp2, SIZE(asca,1),&
         eigvec, SIZE(asca,1),nb,mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1))
#else
    CALL CPP_mult ('L', 'N',m, num2, tmp2, SIZE(asca,1),&
         eigvec, SIZE(asca,1),nb,mpi_comm_rows, mpi_comm_cols, asca, SIZE(asca,1))
#endif

#if defined (CPP_ELPA_201705003)
    CALL elpa_deallocate(elpa_obj)
    CALL elpa_uninit()
#endif
    ! END of ELPA stuff
    CALL BLACS_GRIDEXIT(ictextblacs,ierr)
#if ( !defined (CPP_ELPA_201705003))
    CALL MPI_COMM_FREE(mpi_comm_rows,ierr)
    CALL MPI_COMM_FREE(mpi_comm_cols,ierr)
#endif
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
#ifdef CPP_REAL
    CALL subredist2(m,num2,myrowssca,SUB_COMM,nprow,npcol, myid,ierr,nb,z,asca)
#else
    CALL subredist2(m,num2,myrowssca,SUB_COMM,nprow,npcol, myid,ierr,nb,achi_c=z,asca_c=asca)
#endif
    !

    DEALLOCATE ( asca )
    DEALLOCATE ( bsca )
    DEALLOCATE (eigvec, tmp2, eig2)
 
