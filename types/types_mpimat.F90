!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_mpimat
  USE m_judft
  USE m_types_rcmat
  PRIVATE
  INTEGER,PARAMETER    :: DEFAULT_BLOCKSIZE=64
  INTEGER, PARAMETER   :: dlen_=9

  !<This data-type extends the basic t_mat for distributed matrices.
  !<
  !<It stores the additional mpi_communicator and sets up a blacs grid for the matrix.
  !<This can be used to perform scalapack calls on the matrix with little additional input.
  !<The copy procedure is overwritten from t_mat to enable also redistribution of the matrix.
  
  
  TYPE,EXTENDS(t_mat):: t_mpimat
     INTEGER:: mpi_com                          !> mpi-communiator over which matrix is distributed
     INTEGER:: blacs_desc(dlen_)                !> blacs descriptor
     INTEGER:: blacs_ctext                      !> blacs context
     INTEGER:: global_size1,global_size2        !> this is the size of the full-matrix
     INTEGER:: npcol,nprow                      !> the number of columns/rows in the processor grid
   CONTAINS
     PROCEDURE,PASS   :: copy => mpimat_copy     !<overwriten from t_mat, also performs redistribution
     PROCEDURE,PASS   :: free => mpimat_free     !<overwriten from t_mat, takes care of blacs-grids
     PROCEDURE,PASS   :: init => mpimat_init     !<overwriten from t_mat, also calls alloc in t_mat
     PROCEDURE,PASS   :: add_transpose => mpimat_add_transpose !<overwriten from t_mat
  END TYPE t_mpimat
  
  PUBLIC t_mpimat

CONTAINS

  SUBROUTINE mpimat_add_transpose(mat,mat1)
    CLASS(t_mpimat),INTENT(INOUT) ::mat
    CLASS(t_mat),INTENT(INOUT) ::mat1

    INTEGER:: i,ii,n_size,n_rank

    SELECT TYPE(mat1)
    TYPE IS (t_mpimat)
    
       IF (mat%l_real) THEN
#ifdef CPP_SCALAPACK          

       CALL pdgeadd('t',mat1%global_size1,mat1%global_size2,1.0,mat1%data_r,1,1,mat1%blacs_desc,1.0,mat%data_r,1,1,mat%blacs_desc)
    ELSE
       CALL pzgeadd('t',mat1%global_size1,mat1%global_size2,CMPLX(1.0,0.0),mat1%data_c,1,1,mat1%blacs_desc,CMPLX(1.0,0.0),mat%data_c,1,1,mat%blacs_desc)
#endif
    END IF
    !Now multiply the diagonal of the matrix by 1/2
#ifdef CPP_MPI    
    CALL MPI_COMM_RANK(mat%mpi_com,n_rank,i)
    CALL MPI_COMM_SIZE(mat%mpi_com,n_size,i)
#endif
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
    CLASS(t_mat),INTENT(INOUT)   ::mat1
    INTEGER,INTENT(IN) ::n1,n2
#ifdef CPP_SCALAPACK
    SELECT TYPE(mat1)
    TYPE IS(t_mpimat)
       IF (mat%l_real) THEN
          CALL pdgemr2d(mat1%matsize1,mat1%matsize2,mat1%data_r,1,1,mat1%blacs_desc,mat%data_r,n1,n2,mat%blacs_desc,mat%blacs_ctext)
       ELSE
          CALL pzgemr2d(mat1%matsize1,mat1%matsize2,mat1%data_r,1,1,mat1%blacs_desc,mat%data_r,n1,n2,mat%blacs_desc,mat%blacs_ctext)
       END IF
    CLASS DEFAULT
       CALL judft_error("Wrong datatype in copy")
    END SELECT
#endif    
  END SUBROUTINE mpimat_copy
  
  SUBROUTINE mpimat_free(mat)
    IMPLICIT NONE
    CLASS(t_mpimat),INTENT(INOUT) :: mat
    INTEGER :: ierr
    IF (ALLOCATED(mat%data_r)) DEALLOCATE(mat%data_r)
    IF (ALLOCATED(mat%data_c)) DEALLOCATE(mat%data_c)
#ifdef CPP_SCALAPACK    
    CALL BLACS_GRIDEXIT(mat%blacs_ctext,ierr)
#endif    
  END SUBROUTINE mpimat_free

  !>Initialization of the distributed matrix.
  !!
  !! The argument l_2d controls the kind of distribution used:
  !!  - TRUE: the matrix is a Scalapack BLOCK-CYCLIC distribution
  !!  - FALSE: the matrix is distributed in a one-dimensional column cyclic distribution with blocksize 1
  !! as used in the parallel matrix setup of FLEUR
  SUBROUTINE mpimat_init(mat,l_real,matsize1,matsize2,mpi_subcom,l_2d)
    IMPLICIT NONE
    CLASS(t_mpimat)             :: mat
    INTEGER,INTENT(IN),OPTIONAL :: matsize1,matsize2,mpi_subcom
    LOGICAL,INTENT(IN),OPTIONAL :: l_real,l_2d
    
    INTEGER::nb
    nb=DEFAULT_BLOCKSIZE
    IF (.NOT.(PRESENT(matsize1).AND.PRESENT(matsize2).AND.PRESENT(mpi_subcom).AND.PRESENT(l_real).AND.PRESENT(l_2d)))&
         CALL judft_error("Optional arguments must be present in mpimat_init")
    mat%global_size1=matsize1
    mat%global_size2=matsize2
    mat%mpi_com=mpi_subcom
    CALL priv_create_blacsgrid(mat%mpi_com,l_2d,matsize1,matsize2,nb,&
         mat%blacs_ctext,mat%blacs_desc,&
         mat%matsize1,mat%matsize2,&
         mat%npcol,mat%nprow)
    print *,"mat:",mat%matsize1,mat%matsize2 
    print *,"pe:",mat%npcol,mat%nprow 
    CALL mat%alloc(l_real) !Attention,sizes determined in call to priv_create_blacsgrid
  END SUBROUTINE mpimat_init
    
  SUBROUTINE priv_create_blacsgrid(mpi_subcom,l_2d,m1,m2,nb,ictextblacs,sc_desc,local_size1,local_size2,npcol,nprow)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: mpi_subcom
    INTEGER,INTENT(IN) :: m1,m2
    INTEGER,INTENT(INOUT)::nb
    LOGICAL,INTENT(IN) :: l_2d
    INTEGER,INTENT(OUT):: ictextblacs,sc_desc(:)
    INTEGER,INTENT(OUT):: local_size1,local_size2
    INTEGER,INTENT(OUT):: npcol,nprow
   
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER     :: myrowssca,mycolssca,myrow,mycol
    INTEGER     :: iamblacs,npblacs,np,myid
    INTEGER     :: nprow2,npcol2,myrowblacs,mycolblacs
    INTEGER     :: k,i,j
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
       nb=1
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
    myrowssca=(m1-1)/(nb*nprow)*nb+ MIN(MAX(m1-(m1-1)/(nb*nprow)*nb*nprow-nb*myrow,0),nb)
    !     Number of rows the local process gets in ScaLAPACK distribution
    mycolssca=(m2-1)/(nb*npcol)*nb+ MIN(MAX(m2-(m2-1)/(nb*npcol)*nb*npcol-nb*mycol,0),nb)

    !Get BLACS ranks for all MPI ranks
    CALL BLACS_PINFO(iamblacs,npblacs)  ! iamblacs = local process rank (e.g. myid)
    ! npblacs  = number of available processes
    iblacsnums=-2
    ihelp=-2
    ihelp(myid+1)=iamblacs ! Get the Blacs id corresponding to the MPI id
    !print *,"ALLREDUCE:",mpi_subcom
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
    CALL descinit(sc_desc,m1,m2,nb,nb,0,0,ictextblacs,myrowssca,ierr)
    IF (ierr /=0 ) call judft_error('Creation of BLACS descriptor failed')
    local_size1=myrowssca
    local_size2=mycolssca
#endif
  END SUBROUTINE priv_create_blacsgrid
END MODULE m_types_mpimat
