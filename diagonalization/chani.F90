!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_chani
CONTAINS
  SUBROUTINE chani(M,Neigd, Myid,Np,Sub_comm,mpi_comm, &
       Eig,Num,hamOvlp,zMat)
    !
    !----------------------------------------------------
    !- Parallel eigensystem solver - driver routine; gb99
    !  Uses the SCALAPACK for the actual diagonalization
    !
    ! m ........ actual (=leading) dimension of full a & b matrices
    !            must be problem size, as input a, b  are one-dimensional
    !            and shall be redistributed to two-dimensional matrices
    !            actual (=leading) dimension of eigenvector z(,)
    ! n ........ number of columns of full (sub)matrix ( about n/np)
    ! neigd..... second (=column) dimension of eigenvector matrix
    ! myid ..... id of node (normally irank)
    ! np ....... number of processors (normally isize)
    ! SUB_COMM.. communicator for MPI
    ! a,b   .... packed (sub)matrices, here expanded to non-packed
    ! z,eig .... eigenvectors and values, output
    ! num ...... number of ev's searched (and found) on this node
    !            On input, overall number of ev's searched,
    !            On output, local number of ev's found
    !
    !----------------------------------------------------
    !
!#include"cpp_arch.h"
#include"cpp_double.h"
    USE m_juDFT
    USE m_types
    USE m_subredist1
    USE m_subredist2
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER, INTENT (IN) :: neigd,m
    INTEGER, INTENT (IN) :: SUB_COMM,np,myid,mpi_comm
    INTEGER, INTENT (INOUT) :: num
    REAL,    INTENT   (OUT) :: eig(neigd)
    TYPE(t_hamOvlp),INTENT(INOUT) :: hamOvlp
    TYPE(t_zMat),INTENT(INOUT)    :: zMat
 
    !...  Local variables
    !
    INTEGER nc,ic,ir,n_sym,jsym,num_j,icm,n_bound
    INTEGER i ,j ,k ,l, info, i_plus, i_dim
    INTEGER nprow,npcol,myrowssca, nsq, nprow2, npcol2
    INTEGER myrow,mycol,mycolssca, ierr, me,ierr2
    INTEGER iamblacs,myrowblacs,mycolblacs,npblacs
    INTEGER ictxtblacs,err
    INTEGER, ALLOCATABLE :: iwork(:),iusermap(:,:)
    INTEGER, ALLOCATABLE :: iblacsnums(:),ihelp(:)
    REAL,    ALLOCATABLE :: dwork(:)
    REAL,    ALLOCATABLE :: eigvec_r(:,:)
    COMPLEX, ALLOCATABLE :: eigvec_c(:,:)
    REAL,    ALLOCATABLE :: rwork(:)
    INTEGER              :: lrwork

    !
    !  ScaLAPACK things
    CHARACTER (len=1)    :: uplo
    INTEGER, PARAMETER   :: dlen_=9
    INTEGER              :: desca(dlen_),desceigv(dlen_)
    INTEGER              :: num1,num2,liwork,lwork2,np0,mq0
    INTEGER              :: iceil, numroc, nn, nb
    INTEGER, ALLOCATABLE :: ifail(:), iclustr(:)
    REAL                 :: abstol,orfac=1.E-4,CPP_LAPACK_slamch
    REAL,ALLOCATABLE     :: eig2(:), gap(:)
    REAL,    ALLOCATABLE :: asca_r(:,:), bsca_r(:,:),work2_r(:)
    COMPLEX, ALLOCATABLE :: asca_c(:,:), bsca_c(:,:),work2_c(:)

    EXTERNAL iceil, numroc
    EXTERNAL CPP_LAPACK_slamch, descinit
    EXTERNAL blacs_pinfo, blacs_gridinit
    EXTERNAL MPI_COMM_DUP
    
    !
    ! determine actual number of columns of input matrices A and B
    ! nc is number of columns the local processor will get, must be <=n
    !
    nc = 0
    DO i = myid+1, m, np
       nc = nc + 1
    ENDDO
    !IF (nc.GT.n) THEN
    !   WRITE (6,*) myid,'gets more columns than allowed:'
    !   WRITE (6,*) myid,'will get',nc,' columns, only',n,' allowed'
    !   CALL juDFT_error("chani: nc > n",calledby ="chani")
    !ENDIF
    !
    ! determine block size
    !
    nb = 20
    IF (m.GT.2048)   nb = 30 !2
    IF (m.GT.8*2048) nb = 60 !4
    !
    ! transfer to ScaLAPACK block cyclic 2d distribution with
    ! block-size nb
    !
    ! compute processor grid, as square as possible
    ! If not square with more rows than columns
    nsq = INT(SQRT(REAL(np)))
    gridloop:DO j=nsq,1,-1
       IF ( (np/j) * j == np) THEN
          npcol = np/j
          nprow = j
          EXIT gridloop
       ENDIF
    ENDDO gridloop
    !
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
    myrowssca=(m-1)/(nb*nprow)*nb+&
         MIN(MAX(m-(m-1)/(nb*nprow)*nb*nprow-nb*myrow,0),nb)
    !     Number of rows the local process gets in ScaLAPACK distribution
    mycolssca=(m-1)/(nb*npcol)*nb+ &
         MIN(MAX(m-(m-1)/(nb*npcol)*nb*npcol-nb*mycol,0),nb)
    !     Number of columns the local process gets in ScaLAPACK distribution
    if (hamOvlp%l_real) THEN
       ALLOCATE ( asca_r(myrowssca,mycolssca), stat=err )
    else
       ALLOCATE ( asca_c(myrowssca,mycolssca), stat=err )
    endif
    IF (err.NE.0) THEN
       WRITE (*,*) 'In chani an error occured during the allocation of'
       WRITE (*,*) 'asca: ',err,' size: ',myrowssca*mycolssca
       CALL juDFT_error("allocte in chani",calledby ="chani")
    ENDIF
    if (hamOvlp%l_real) THEN
       asca_r=0.0
       ALLOCATE ( bsca_r(myrowssca,mycolssca), stat=err  )
    else
       asca_c=0.0
       ALLOCATE ( bsca_c(myrowssca,mycolssca), stat=err  )
    endif
    IF (err.NE.0) THEN
       WRITE (*,*) 'In chani an error occured during the allocation of'
       WRITE (*,*) 'bsca :',err
       CALL juDFT_error("allocate in chani",calledby ="chani")
    ENDIF
    if (hamOvlp%l_real) THEN
       bsca_r=0.0
       CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb ,achi_r=hamovlp%a_r,asca_r=asca_r)
       CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb ,achi_r=hamovlp%b_r,asca_r=bsca_r)
    else
       bsca_c=0.0
       CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb ,achi_c=hamovlp%a_c,asca_c=asca_c)
       CALL subredist1(m,myrowssca,SUB_COMM, nprow, npcol, myid, ierr, nb ,achi_c=hamovlp%b_c,asca_c=bsca_c)
    end if
    CALL BLACS_PINFO(iamblacs,npblacs)  ! iamblacs = local process rank (e.g. myid)
    ! npblacs  = number of available processes
    IF (npblacs /= np) THEN
       !gb         WRITE (6,*) 'Number of procs differ between',
       !gb     &        ' BLACS and MPI'
       IF (np > npblacs) THEN
          WRITE(6,*) 'Have to take SUB_COMM for BLACS grid!' ! this should not happen
       ENDIF
    ENDIF
    ALLOCATE( iblacsnums(np) )
    iblacsnums=-2
    ALLOCATE (ihelp(np))
    ihelp=-2
    ihelp(myid+1)=iamblacs ! Get the Blacs id corresponding to the MPI id
    CALL MPI_ALLREDUCE(ihelp, iblacsnums, np,MPI_INTEGER,MPI_MAX,SUB_COMM,ierr)
    IF (ierr.NE.0) THEN
       WRITE (6,*) 'Error in allreduce for BLACS nums'
       CALL juDFT_error('Error in allreduce for BLACS nums',calledby ='chani')
    ENDIF
    !     iblacsnums(i) is the BLACS-process number of MPI-process i-1
    ALLOCATE ( iusermap(nprow,npcol) ) ! Usermap for BLACS_GRIDMAP
    k = 1
    DO i = 1, nprow
       DO j = 1, npcol
          iusermap(i,j) = iblacsnums(k)
          k = k + 1
       ENDDO
    ENDDO
    !gb      WRITE(*,*) np,iblacsnums
    !
    !     Take a copy of SUB_COMM as system context for BLACS_GRIDMAP
    !     We cannot take SUB_COMM itself as it will be overwritten with the
    !     BLACS context for the grid comstructed 
    !
    !      CALL MPI_COMM_DUP(SUB_COMM,ictxtblacs,ierr)
    !I do not understand this, perhaps we should use MPI_COMM_WRLD???
    !      CALL MPI_COMM_DUP(MPI_COMM_WRLD,ictxtblacs,ierr)
    CALL MPI_COMM_DUP(MPI_COMM,ictxtblacs,ierr)
    IF (ierr /=0 ) THEN
       WRITE(6,*) 'MPI_COMM_DUP failed'
       CALL juDFT_error('MPI_COMM_DUP failed',calledby='chani')
    ENDIF
    CALL BLACS_GRIDMAP(ictxtblacs,iusermap,nprow,nprow,npcol)
    !     Now control, whether the BLACS grid is the one we wanted
    CALL BLACS_GRIDINFO(ictxtblacs, nprow2,npcol2,myrowblacs,mycolblacs)
    IF (nprow2 /= nprow) THEN
       WRITE(6,*) 'Wrong number of rows in BLACS grid'
       WRITE(6,*) 'nprow=',nprow,' nprow2=',nprow2
       ierr = -nprow2
       CALL juDFT_error('Wrong number of rows in BLACS grid',calledby= 'chani')
    ENDIF
    IF (npcol2 /= npcol) THEN
       WRITE(6,*) 'Wrong number of columns in BLACS grid'
       WRITE(6,*) 'npcol=',npcol,' npcol2=',npcol2
       ierr = -npcol2
       CALL juDFT_error('Wrong number of columns in BLACS grid', calledby ='chani')

    ENDIF
    IF (myrowblacs.NE.myrow) THEN
       WRITE(6,*) 'My row in BLACS is',myrowblacs
       WRITE(6,*) 'It is not',myrow
       !CALL CPP_flush(6)
    ENDIF
    IF (mycolblacs.NE.mycol) THEN
       WRITE(6,*) 'My column in BLACS is',mycolblacs
       WRITE(6,*) 'It is not',mycol
       !CALL CPP_flush(6)
    ENDIF
    !
    CALL descinit(DESCA,m,m,nb,nb,0,0,ictxtblacs,myrowssca,ierr)
    CALL descinit(DESCeigv,m,m,nb,nb,0,0,ictxtblacs, myrowssca,ierr)

    abstol=2.0*CPP_LAPACK_slamch('S') ! PDLAMCH gave an error on ZAMpano
    ALLOCATE ( eig2(m), stat=err ) ! The eigenvalue array for ScaLAPACK
    IF (err.NE.0) THEN
       WRITE (*,*) 'In chani an error occured during the allocation of'
       WRITE (*,*) 'eig2 :',err
       CALL juDFT_error('Failed to allocated "eig2"', calledby ='chani')
    ENDIF
    !      write(*,*) 'c :',myrowssca,mycolssca,desceigv
    if (hamovlp%l_real) THEN
       ALLOCATE ( eigvec_r(myrowssca,mycolssca), stat=err ) ! Eigenvectors for ScaLAPACK
    else
       ALLOCATE ( eigvec_c(myrowssca,mycolssca), stat=err ) ! Eigenvectors for ScaLAPACK
    end if
    IF (err.NE.0) THEN
       WRITE (*,*) 'In chani an error occured during the allocation of'
       WRITE (*,*) 'eigvec :',err
       CALL juDFT_error('Failed to allocated "eigvec"', calledby ='chani')
    ENDIF
    !
    nn=MAX(MAX(m,nb),2)
    np0=numroc(nn,nb,0,0,nprow)
    mq0=numroc(MAX(MAX(neigd,nb),2),nb,0,0,npcol)
    if (hamovlp%l_real) THEN
       lwork2=5*m+MAX(5*nn,np0*mq0+2*nb*nb)+ iceil(neigd,nprow*npcol)*nn
       ALLOCATE ( work2_r(lwork2+10*m), stat=err ) ! Allocate more in case of clusters
    else
       lwork2=m+MAX(nb*(np0+1),3)
       ALLOCATE ( work2_c(lwork2), stat=err )
    endif
    IF (err.NE.0) THEN 
       WRITE (*,*) 'work2  :',err,lwork2
       CALL juDFT_error('Failed to allocated "work2"', calledby ='chani')
    ENDIF

    liwork=6*MAX(MAX(m,nprow*npcol+1),4)
    ALLOCATE ( iwork(liwork), stat=err )
    IF (err.NE.0) THEN
       WRITE (*,*) 'iwork  :',err,liwork
       CALL juDFT_error('Failed to allocated "iwork"', calledby ='chani')
    ENDIF
    ALLOCATE ( ifail(m), stat=err )
    IF (err.NE.0) THEN
       WRITE (*,*) 'ifail  :',err,m
       CALL juDFT_error('Failed to allocated "ifail"', calledby ='chani')
    ENDIF
    ALLOCATE ( iclustr(2*nprow*npcol), stat=err )
    IF (err.NE.0) THEN
       WRITE (*,*) 'iclustr:',err,2*nprow*npcol
       CALL juDFT_error('Failed to allocated "iclustr"', calledby ='chani')
    ENDIF
    ALLOCATE ( gap(nprow*npcol), stat=err )
    IF (err.NE.0) THEN
       WRITE (*,*) 'gap    :',err,nprow*npcol
       CALL juDFT_error('Failed to allocated "gap"', calledby ='chani')
    ENDIF
    !
    !     Compute size of workspace
    !
    if (hamovlp%l_real) THEN
    uplo='U'
    CALL CPP_LAPACK_pdsygvx(1,'V','I','U',m,asca_r,1,1,desca,bsca_r,1,1, desca,&
         0.0,1.0,1,num,abstol,num1,num2,eig2,orfac,eigvec_r,1,1,&
         desceigv,work2_r,-1,iwork,-1,ifail,iclustr, gap,ierr)
    IF ( work2_r(1).GT.lwork2) THEN
       lwork2 = work2_r(1)
       DEALLOCATE (work2_r)
       ALLOCATE ( work2_r(lwork2+20*m), stat=err ) ! Allocate even more in case of clusters
       IF (err.NE.0) THEN
          WRITE (*,*) 'work2  :',err,lwork2
          CALL juDFT_error('Failed to allocated "work2"', calledby ='chani')
       ENDIF
    ENDIF
else
    lrwork=4*m+MAX(5*nn,np0*mq0)+ iceil(neigd,nprow*npcol)*nn
    ! Allocate more in case of clusters
    ALLOCATE(rwork(lrwork+10*m), stat=ierr)
    IF (err /= 0) THEN
       WRITE (*,*) 'ERROR: chani.F: Allocating rwork failed'
       CALL juDFT_error('Failed to allocated "rwork"', calledby ='chani')
    ENDIF

    CALL CPP_LAPACK_pzhegvx(1,'V','I','U',m,asca_c,1,1,desca,bsca_c,1,1, desca,&
         0.0,1.0,1,num,abstol,num1,num2,eig2,orfac,eigvec_c,1,1,&
         desceigv,work2_c,-1,rwork,-1,iwork,-1,ifail,iclustr,&
         gap,ierr)
    IF (ABS(work2_c(1)).GT.lwork2) THEN
       lwork2=work2_c(1)
       DEALLOCATE (work2_c)
       ALLOCATE (work2_c(lwork2), stat=err)
       IF (err /= 0) THEN
          WRITE (*,*) 'ERROR: chani.F: Allocating rwork failed:',lwork2
          CALL juDFT_error('Failed to allocated "work2"', calledby ='chani')
       ENDIF
    ENDIF
    IF (rwork(1).GT.lrwork) THEN
       lrwork=rwork(1)
       DEALLOCATE(rwork)
       ! Allocate even more in case of clusters
       ALLOCATE (rwork(lrwork+20*m), stat=err)
       IF (err /= 0) THEN
          WRITE (*,*) 'ERROR: chani.F: Allocating rwork failed: ', lrwork+20*m
          CALL juDFT_error('Failed to allocated "rwork"', calledby ='chani')
       ENDIF
    ENDIF
endif
    IF (iwork(1) .GT. liwork) THEN
       liwork = iwork(1)
       DEALLOCATE (iwork)
       ALLOCATE (iwork(liwork), stat=err)
       IF (err /= 0) THEN
          WRITE (*,*) 'ERROR: chani.F: Allocating iwork failed: ',liwork
          CALL juDFT_error('Failed to allocated "iwork"', calledby ='chani')
       ENDIF
    ENDIF
    !
    !     Now solve generalized eigenvalue problem
    !
if (hamovlp%l_real) THEN
    CALL CPP_LAPACK_pdsygvx(1,'V','I','U',m,asca_r,1,1,desca,bsca_r,1,1, desca,&
         1.0,1.0,1,num,abstol,num1,num2,eig2,orfac,eigvec_r,1,1,&
         desceigv,work2_r,lwork2,iwork,liwork,ifail,iclustr,&
         gap,ierr)
else
    CALL CPP_LAPACK_pzhegvx(1,'V','I','U',m,asca_c,1,1,desca,bsca_c,1,1, desca,&
         1.0,1.0,1,num,abstol,num1,num2,eig2,orfac,eigvec_c,1,1,&
         desceigv,work2_c,lwork2,rwork,lrwork,iwork,liwork,&
         ifail,iclustr,gap,ierr)
    DEALLOCATE(rwork)
endif
    IF (ierr .NE. 0) THEN
       IF (ierr /= 2) WRITE (6,*) myid,' error in pzhegvx/pdsygvx, ierr=',ierr
       IF (ierr <= 0) WRITE (6,*) myid,' illegal input argument'
       IF (MOD(ierr,2) /= 0) THEN
          WRITE (6,*) myid,'some eigenvectors failed to converge'
          eigs: DO i = 1, neigd
             IF (ifail(i) /= 0) THEN
                WRITE (6,*) myid,' eigenvector',ifail(i), 'failed to converge'
             ELSE
                EXIT eigs
             ENDIF
          ENDDO eigs
          !CALL CPP_flush(6)
       ENDIF
       IF (MOD(ierr/4,2).NE.0) THEN
          WRITE(6,*) myid,' only',num2,' eigenvectors converged'
          !CALL CPP_flush(6)
       ENDIF
       IF (MOD(ierr/8,2).NE.0) THEN
          WRITE(6,*) myid,' PDSTEBZ failed to compute eigenvalues'
          CALL judft_error("SCALAPACK failed to solve eigenvalue problem",calledby="chani.F90")
       ENDIF
       IF (MOD(ierr/16,2).NE.0) THEN
          WRITE(6,*) myid,' B was not positive definite, Cholesky failed at',ifail(1)
          CALL judft_error("SCALAPACK failed: B was not positive definite",calledby="chani.F90")
       ENDIF
    ENDIF
    IF (num2 < num1) THEN
       IF (myid ==0) THEN
          WRITE(6,*) 'Not all eigenvalues wanted are found'
          WRITE(6,*) 'number of eigenvalues/vectors wanted',num1
          WRITE(6,*) 'number of eigenvalues/vectors found',num2
          !CALL CPP_flush(6)
       ENDIF
    ENDIF
    ierr=0
    CALL BLACS_GRIDEXIT(ictxtblacs,ierr)
    IF (ierr /= 0 ) THEN
       WRITE (6,*) 'blacs_gridexit failed, ierr=',ierr
    ENDIF
    DEALLOCATE(iwork)
    DEALLOCATE(ifail)
    DEALLOCATE(iclustr)
    DEALLOCATE(gap)
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
    DEALLOCATE(eig2)
    !
    !     Redistribute eigvec from ScaLAPACK distribution to each process
    !     having all eigenvectors corresponding to his eigenvalues as above
    !
    if (hamovlp%l_real) THEN
       CALL subredist2(m,num2,myrowssca,SUB_COMM,nprow,npcol, myid,ierr,nb,zmat%z_r,eigvec_r)
       ELSE
       CALL subredist2(m,num2,myrowssca,SUB_COMM,nprow,npcol, myid,ierr,nb,achi_c=zmat%z_c,asca_c=eigvec_c)
    end if
    !
    !DEALLOCATE ( eigvec)
    DEALLOCATE ( iblacsnums )
    DEALLOCATE ( ihelp )
    DEALLOCATE ( iusermap )

  END SUBROUTINE chani
end MODULE m_chani
