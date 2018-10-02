!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_mpimat
  USE m_judft
  USE m_types_rcmat
  IMPLICIT NONE
  PRIVATE
  INTEGER,PARAMETER    :: DEFAULT_BLOCKSIZE=64
  INTEGER, PARAMETER   :: dlen_=9

  !<This data-type extends the basic t_mat for distributed matrices.
  !<
  !<It stores the additional mpi_communicator and sets up a blacs grid for the matrix.
  !<This can be used to perform scalapack calls on the matrix with little additional input.
  !<The copy procedure is overwritten from t_mat to enable also redistribution of the matrix.
  TYPE  t_blacsdata
     INTEGER:: no_use
     INTEGER:: mpi_com                          !> mpi-communiator over which matrix is distributed
     INTEGER:: blacs_desc(dlen_)                !> blacs descriptor
     !> 1: =1
     !> 2: context
     !> 3,4: global matrix size
     !> 5,6: block sizes
     !> 7,8: row/colum of grid for first row/colum of matrix
     !> 9: leading dimension of local matrix
     INTEGER:: npcol,nprow                     !> the number of columns/rows in the processor grid
  END TYPE t_blacsdata
  
  
  TYPE,EXTENDS(t_mat):: t_mpimat
     INTEGER                   :: global_size1,global_size2        !> this is the size of the full-matrix
     TYPE(t_blacsdata),POINTER :: blacsdata
   CONTAINS
     PROCEDURE,PASS   :: copy => mpimat_copy     !<overwriten from t_mat, also performs redistribution
     PROCEDURE,PASS   :: move => mpimat_move     !<overwriten from t_mat, also performs redistribution
     PROCEDURE,PASS   :: free => mpimat_free     !<overwriten from t_mat, takes care of blacs-grids
     PROCEDURE,PASS   :: init_details => mpimat_init
     PROCEDURE,PASS   :: init_template =>mpimat_init_template     !<overwriten from t_mat, also calls alloc in t_mat
     PROCEDURE,PASS   :: add_transpose => mpimat_add_transpose !<overwriten from t_mat
     PROCEDURE,PASS   :: generate_full_matrix    ! construct full matrix if only upper triangle of hermitian matrix is given
     PROCEDURE,PASS   :: print_matrix
     PROCEDURE,PASS   :: from_non_dist
     FINAL :: finalize
  END TYPE t_mpimat
  
  PUBLIC t_mpimat

CONTAINS

  SUBROUTINE print_matrix(mat,fileno)
    CLASS(t_mpimat),INTENT(INOUT) ::mat
    INTEGER:: fileno

#ifdef CPP_SCALAPACK
    INCLUDE 'mpif.h'
    INTEGER,EXTERNAL:: indxl2g
    CHARACTER(len=10)::filename
    INTEGER :: irank,isize,i,j,npr,npc,r,c,tmp,err,status(MPI_STATUS_SIZE) 

    CALL MPI_COMM_RANK(mat%blacsdata%mpi_com,irank,err)
    CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com,isize,err)

    tmp=0

    IF (irank>0) CALL MPI_RECV(tmp,1,MPI_INTEGER,irank-1,0,mat%blacsdata%mpi_com,status,err) !lock
    WRITE(filename,"(a,i0)") "out.",fileno
    OPEN(fileno,file=filename,access='append')
    
    CALL blacs_gridinfo(mat%blacsdata%blacs_desc(2),npr,npc,r,c)
    DO i=1,mat%matsize1
       DO j=1,mat%matsize2
          IF (mat%l_real) THEN
             WRITE(fileno,"(5(i0,1x),2(f10.5,1x))") irank,i,j,indxl2g(i,mat%blacsdata%blacs_desc(5),r,0,npr),&
                  indxl2g(j,mat%blacsdata%blacs_desc(6),c,0,npc),mat%data_r(i,j)
          ELSE
             WRITE(fileno,"(5(i0,1x),2(f10.5,1x))") irank,i,j,indxl2g(i,mat%blacsdata%blacs_desc(5),r,0,npr),&
                  indxl2g(j,mat%blacsdata%blacs_desc(6),c,0,npc),mat%data_c(i,j)
          END IF
       ENDDO
    ENDDO
    CLOSE(fileno)
    IF (irank+1<isize) CALL MPI_SEND(tmp,1,MPI_INTEGER,irank+1,0,mat%blacsdata%mpi_com,err)
    
#endif    
  END SUBROUTINE print_matrix

  SUBROUTINE generate_full_matrix(mat)
    CLASS(t_mpimat),INTENT(INOUT) ::mat
    
    INTEGER :: i,j,i_glob,j_glob,npcol,nprow,myid,err,myrow,mycol,np
    COMPLEX,ALLOCATABLE:: tmp_c(:,:)
    REAL,ALLOCATABLE   :: tmp_r(:,:)
#ifdef CPP_SCALAPACK
    INCLUDE 'mpif.h'
    INTEGER, EXTERNAL    :: numroc, indxl2g  !SCALAPACK functions

    CALL MPI_COMM_RANK(mat%blacsdata%mpi_com,myid,err)
    CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com,np,err)
 
    CALL blacs_gridinfo(mat%blacsdata%blacs_desc(2),nprow,npcol,myrow,mycol)
    !Set lower part of matrix to zero
 
    DO i=1,mat%matsize1
       DO j=1,mat%matsize2
          ! Get global column corresponding to i and number of local rows up to
          ! and including the diagonal, these are unchanged in A
          i_glob = indxl2g(i,     mat%blacsdata%blacs_desc(5), myrow, 0, nprow)
          j_glob = indxl2g(j,     mat%blacsdata%blacs_desc(6), mycol, 0, npcol)

          IF (i_glob>j_glob) THEN
             IF (mat%l_real) THEN
                mat%data_r(i,j) = 0.0
             ELSE
                mat%data_c(i,j) = 0.0
             ENDIF
          ENDIF
          IF (i_glob==j_glob) THEN
             IF (mat%l_real) THEN
                mat%data_r(i,j) =  mat%data_r(i,j)/2.0
             ELSE
                mat%data_c(i,j) =  mat%data_c(i,j)/2.0
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    IF (mat%l_real) THEN
       ALLOCATE(tmp_r(mat%matsize1,mat%matsize2))
       tmp_r=mat%data_r
    ELSE
       ALLOCATE(tmp_c(mat%matsize1,mat%matsize2))
       tmp_c=mat%data_c
    END IF
    CALL MPI_BARRIER(mat%blacsdata%mpi_com,i)
 IF (mat%l_real) THEN
#ifdef CPP_SCALAPACK          

       CALL pdgeadd('t',mat%global_size1,mat%global_size2,1.0,tmp_r,1,1,mat%blacsdata%blacs_desc,1.0,mat%data_r,1,1,mat%blacsdata%blacs_desc)
    ELSE
       CALL pzgeadd('c',mat%global_size1,mat%global_size2,CMPLX(1.0,0.0),tmp_c,1,1,mat%blacsdata%blacs_desc,CMPLX(1.0,0.0),mat%data_c,1,1,mat%blacsdata%blacs_desc)
#endif
    END IF


#endif
  END SUBROUTINE generate_full_matrix


  SUBROUTINE mpimat_add_transpose(mat,mat1)
    CLASS(t_mpimat),INTENT(INOUT) ::mat
    CLASS(t_mat),INTENT(INOUT) ::mat1

    INTEGER:: i,ii,n_size,n_rank

    SELECT TYPE(mat1)
    TYPE IS (t_mpimat)
#ifdef CPP_MPI    
    CALL MPI_COMM_RANK(mat%blacsdata%mpi_com,n_rank,i)
    CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com,n_size,i)
#endif
    !Set lower part of matrix to zero...
       ii=0
       DO i=n_rank+1,MIN(mat%global_size1,mat%global_size2),n_size
          ii=ii+1
          IF (mat%l_real) THEN
             mat%data_r(i+1:,ii)=0.0
             mat1%data_r(i+1:,ii)=0.0
          ELSE
             mat%data_c(i+1:,ii)=0.0
             mat1%data_c(i+1:,ii)=0.0
          ENDIF
       ENDDO
       IF (mat%l_real) THEN
#ifdef CPP_SCALAPACK          

       CALL pdgeadd('t',mat1%global_size1,mat1%global_size2,1.0,mat1%data_r,1,1,mat1%blacsdata%blacs_desc,1.0,mat%data_r,1,1,mat%blacsdata%blacs_desc)
    ELSE
       CALL pzgeadd('c',mat1%global_size1,mat1%global_size2,CMPLX(1.0,0.0),mat1%data_c,1,1,mat1%blacsdata%blacs_desc,CMPLX(1.0,0.0),mat%data_c,1,1,mat1%blacsdata%blacs_desc)
#endif
    END IF
    !Now multiply the diagonal of the matrix by 1/2

    ii=0
    DO i=n_rank+1,MIN(mat%global_size1,mat%global_size2),n_size
       ii=ii+1
       IF (mat%l_real) THEN
          mat%data_r(i,ii)=mat%data_r(i,ii)/2
       ELSE
          mat%data_c(i,ii)=mat%data_c(i,ii)/2
       END IF
    ENDDO
    CLASS default
       CALL judft_error("Inconsistent types in t_mpimat_add_transpose")
    END SELECT
    
  END SUBROUTINE mpimat_add_transpose

  SUBROUTINE mpimat_copy(mat,mat1,n1,n2)
    IMPLICIT NONE
    CLASS(t_mpimat),INTENT(INOUT)::mat
    CLASS(t_mat),INTENT(IN)      ::mat1
    INTEGER,INTENT(IN) ::n1,n2
#ifdef CPP_SCALAPACK
    SELECT TYPE(mat1)
    TYPE IS(t_mpimat)
       IF (mat%l_real) THEN
          CALL pdgemr2d(Mat1%global_size1,mat1%global_size2,mat1%data_r,1,1,mat1%blacsdata%blacs_desc,mat%data_r,n1,n2,mat%blacsdata%blacs_desc,mat%blacsdata%blacs_desc(2))
       ELSE
          CALL pzgemr2d(mat1%global_size1,mat1%global_size2,mat1%data_c,1,1,mat1%blacsdata%blacs_desc,mat%data_c,n1,n2,mat%blacsdata%blacs_desc,mat%blacsdata%blacs_desc(2))
       END IF
    CLASS DEFAULT
       CALL judft_error("Wrong datatype in copy")
    END SELECT
#endif    
  END SUBROUTINE mpimat_copy

  SUBROUTINE from_non_dist(mat,mat1)
    IMPLICIT NONE
    CLASS(t_mpimat),INTENT(INOUT)::mat
    TYPE(t_mat),INTENT(IN)       ::mat1

    INTEGER:: blacs_desc(9),irank,ierr,umap(1,1),np
#ifdef CPP_SCALAPACK
    blacs_desc=(/1,-1,mat1%matsize1,mat1%matsize2,mat1%matsize1,mat1%matsize2,0,0,mat1%matsize1/)

    CALL MPI_COMM_RANK(mat%blacsdata%mpi_com,irank,ierr)
    umap(1,1)=0
    CALL BLACS_GET(mat%blacsdata%blacs_desc(2),10,blacs_desc(2))
    CALL BLACS_GRIDMAP(blacs_desc(2),umap,1,1,1)
    IF (mat%l_real) THEN
       CALL pdgemr2d(Mat1%matsize1,mat1%matsize2,mat1%data_r,1,1,blacs_desc,mat%data_r,1,1,mat%blacsdata%blacs_desc,mat%blacsdata%blacs_desc(2))
    ELSE
       CALL pzgemr2d(mat1%matsize1,mat1%matsize2,mat1%data_c,1,1,blacs_desc,mat%data_c,1,1,mat%blacsdata%blacs_desc,mat%blacsdata%blacs_desc(2))
    END IF
#endif    
  END SUBROUTINE from_non_dist
    


  SUBROUTINE mpimat_move(mat,mat1)
    IMPLICIT NONE
    CLASS(t_mpimat),INTENT(INOUT)::mat
    CLASS(t_mat),INTENT(INOUT)   ::mat1
    CALL mat%copy(mat1,1,1)
  END SUBROUTINE mpimat_move

  SUBROUTINE finalize(mat)
    IMPLICIT NONE
    TYPE(t_mpimat),INTENT(INOUT) :: mat
    CALL mpimat_free(mat)
  END SUBROUTINE finalize

  SUBROUTINE mpimat_free(mat)
    IMPLICIT NONE
    CLASS(t_mpimat),INTENT(INOUT) :: mat
    INTEGER :: ierr
    IF (ALLOCATED(mat%data_r)) DEALLOCATE(mat%data_r)
    IF (ALLOCATED(mat%data_c)) DEALLOCATE(mat%data_c)
    IF (ASSOCIATED(mat%blacsdata)) THEN
       IF (mat%blacsdata%no_use>1) THEN
          mat%blacsdata%no_use=mat%blacsdata%no_use-1
          mat%blacsdata=>null()
       ELSE
#ifdef CPP_SCALAPACK    
          CALL BLACS_GRIDEXIT(mat%blacsdata%blacs_desc(2),ierr)
          DEALLOCATE(mat%blacsdata)
#endif
       END IF
    ENDIF
  END SUBROUTINE mpimat_free

  !>Initialization of the distributed matrix.
  !!
  !! The argument l_2d controls the kind of distribution used:
  !!  - TRUE: the matrix is a Scalapack BLOCK-CYCLIC distribution
  !!  - FALSE: the matrix is distributed in a one-dimensional column cyclic distribution with blocksize 1
  !! as used in the parallel matrix setup of FLEUR
  SUBROUTINE mpimat_init(mat,l_real,matsize1,matsize2,mpi_subcom,l_2d,nb_x,nb_y)
    IMPLICIT NONE
    CLASS(t_mpimat)             :: mat
    INTEGER,INTENT(IN),OPTIONAL :: matsize1,matsize2,mpi_subcom
    LOGICAL,INTENT(IN),OPTIONAL :: l_real,l_2d
    INTEGER,INTENT(IN),OPTIONAL :: nb_y,nb_x
#ifdef CPP_SCALAPACK    
    INTEGER::nbx,nby,irank,ierr
    include 'mpif.h'
    nbx=DEFAULT_BLOCKSIZE; nby=DEFAULT_BLOCKSIZE
    IF (PRESENT(nb_x)) nbx=nb_x
    IF (PRESENT(nb_y)) nby=nb_y
    IF (.NOT.(PRESENT(matsize1).AND.PRESENT(matsize2).AND.PRESENT(mpi_subcom).AND.PRESENT(l_real).AND.PRESENT(l_2d)))&
         CALL judft_error("Optional arguments must be present in mpimat_init")
    mat%global_size1=matsize1
    mat%global_size2=matsize2
    ALLOCATE(mat%blacsdata)
    mat%blacsdata%no_use=1
    mat%blacsdata%mpi_com=mpi_subcom
    CALL priv_create_blacsgrid(mat%blacsdata%mpi_com,l_2d,matsize1,matsize2,nbx,nby,&
         mat%blacsdata%blacs_desc,&
         mat%matsize1,mat%matsize2,&
         mat%blacsdata%npcol,mat%blacsdata%nprow)
    CALL mat%alloc(l_real) !Attention,sizes determined in call to priv_create_blacsgrid
    !check if this matrix is actually distributed over MPI_COMM_SELF
    IF (mpi_subcom==MPI_COMM_SELF) THEN
       CALL MPI_COMM_RANK(mpi_subcom,irank,ierr)
       IF (irank>0) mat%blacsdata%blacs_desc(2)=-1
    END IF
#endif    
  END SUBROUTINE mpimat_init

  SUBROUTINE mpimat_init_template(mat,templ)
    IMPLICIT NONE
    CLASS(t_mpimat),INTENT(INOUT)  :: mat
    CLASS(t_mat),INTENT(IN)        :: templ

    SELECT TYPE(templ)
    TYPE IS (t_mpimat)
       mat%l_real=templ%l_real
       mat%matsize1=templ%matsize1
       mat%matsize2=templ%matsize2
       mat%global_size1=templ%global_size1
       mat%global_size2=templ%global_size2
       mat%blacsdata=>templ%blacsdata
       mat%blacsdata%no_use=mat%blacsdata%no_use+1
       CALL mat%alloc()
      
       CLASS default
          CALL judft_error("Mixed initialization in t_mpimat not possible(BUG)")
    END SELECT
  END SUBROUTINE mpimat_init_template

    
  SUBROUTINE priv_create_blacsgrid(mpi_subcom,l_2d,m1,m2,nbc,nbr,sc_desc,local_size1,local_size2,nprow,npcol)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: mpi_subcom
    INTEGER,INTENT(IN) :: m1,m2
    INTEGER,INTENT(INOUT)::nbc,nbr
    LOGICAL,INTENT(IN) :: l_2d
    INTEGER,INTENT(OUT):: sc_desc(:)
    INTEGER,INTENT(OUT):: local_size1,local_size2
    INTEGER,INTENT(OUT):: npcol,nprow
   
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER     :: myrowssca,mycolssca,myrow,mycol
    INTEGER     :: iamblacs,npblacs,np,myid
    INTEGER     :: nprow2,npcol2,myrowblacs,mycolblacs
    INTEGER     :: k,i,j,ictextblacs
    INTEGER     :: ierr

    INTEGER,ALLOCATABLE :: iblacsnums(:),ihelp(:),iusermap(:,:)

    EXTERNAL descinit, blacs_get
    EXTERNAL blacs_pinfo, blacs_gridinit
 
    !Determine rank and no of processors
    CALL MPI_COMM_RANK(mpi_subcom,myid,ierr)
    CALL MPI_COMM_SIZE(mpi_subcom,np,ierr)


    ! compute processor grid, as square as possible
    ! If not square with more rows than columns
    IF (l_2d) THEN
       distloop: DO j=INT(SQRT(REAL(np))),1,-1
          IF ( (np/j) * j == np) THEN
             npcol = np/j
             nprow = j
             EXIT distloop
          ENDIF
       ENDDO distloop
    ELSE
       nbc=1
       nbr=MAX(m1,m2)
       npcol=np
       nprow=1
    ENDIF
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
    myrowssca=(m1-1)/(nbr*nprow)*nbr+ MIN(MAX(m1-(m1-1)/(nbr*nprow)*nbr*nprow-nbr*myrow,0),nbr)
    !     Number of rows the local process gets in ScaLAPACK distribution
    mycolssca=(m2-1)/(nbc*npcol)*nbc+ MIN(MAX(m2-(m2-1)/(nbc*npcol)*nbc*npcol-nbc*mycol,0),nbc)
  

    !Get BLACS ranks for all MPI ranks
    CALL BLACS_PINFO(iamblacs,npblacs)  ! iamblacs = local process rank (e.g. myid)
    ! npblacs  = number of available processes
    iblacsnums=-2
    ihelp=-2
    ihelp(myid+1)=iamblacs ! Get the Blacs id corresponding to the MPI id
    CALL MPI_ALLREDUCE(ihelp, iblacsnums, np,MPI_INTEGER,MPI_MAX,mpi_subcom,ierr)
    IF (ierr.NE.0) STOP 'Error in allreduce for BLACS nums' 

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
       call judft_error('Wrong number of rows in BLACS grid')
    ENDIF
    IF (npcol2 /= npcol) THEN
       WRITE(6,*) 'Wrong number of columns in BLACS grid'
       WRITE(6,*) 'npcol=',npcol,' npcol2=',npcol2
       call judft_error('Wrong number of columns in BLACS grid')

    ENDIF

    !Create the descriptors
    CALL descinit(sc_desc,m1,m2,nbr,nbc,0,0,ictextblacs,myrowssca,ierr)
    IF (ierr /=0 ) call judft_error('Creation of BLACS descriptor failed')
    local_size1=myrowssca
    local_size2=mycolssca
#endif
  END SUBROUTINE priv_create_blacsgrid
END MODULE m_types_mpimat
