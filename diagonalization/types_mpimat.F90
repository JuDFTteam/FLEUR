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
  
  TYPE,EXTENDS(t_mat):: t_mpimat
     INTEGER:: mpi_com
     INTEGER:: blacs_desc(dlen_)
     INTEGER:: blacs_ctext
     INTEGER:: global_size1,global_size2
     INTEGER:: npcol,nprow

   CONTAINS
     PROCEDURE,PASS   :: init =>mpimat_init
     PROCEDURE,NOPASS :: create_blacsgrid => mpimat_create_blacsgrid
#ifdef CPP_F03
     FINAL            :: free => mpimat_free
#else
     PROCEDURE        :: free => mpimat_free
#endif     
  END TYPE t_mpimat
  
  PUBLIC t_mpimat

CONTAINS

  SUBROUTINE mpimat_free(mat)
    IMPLICIT NONE
    CLASS(t_mpimat) :: mat
    INTEGER :: ierr
    IF (ALLOCATED(mat%data_r)) DEALLOCATE(mat%data_r)
    IF (ALLOCATED(mat%data_c)) DEALLOCATE(mat%data_c)
#ifdef CPP_SCALAPACK    
    CALL BLACS_GRIDEXIT(mat%blacs_ctext,ierr)
#endif    
  END SUBROUTINE mpimat_free
  
  SUBROUTINE mpimat_init(mat,mpi_subcom,m1,m2,l_real,l_2d)
    IMPLICIT NONE
    CLASS(t_mpimat),INTENT(OUT):: mat
    INTEGER,INTENT(IN)         :: m1,m2,mpi_subcom
    LOGICAL,INTENT(IN)         :: l_real,l_2d
    mat%global_size1=m1
    mat%global_size2=m2
    mat%mpi_com=mpi_subcom
    CALL mpimat_create_blacsgrid(mat%mpi_com,l_2d,m1,m2,DEFAULT_BLOCKSIZE,&
         mat%blacs_ctext,mat%blacs_desc,&
         mat%matsize1,mat%matsize2,&
         mat%npcol,mat%nprow)
    CALL mat%alloc(l_real,m1,m2)
    
  END SUBROUTINE mpimat_init

  SUBROUTINE mpimat_create_blacsgrid(mpi_subcom,l_2d,m1,m2,nb,ictextblacs,sc_desc,local_size1,local_size2,npcol,nprow)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: mpi_subcom
    INTEGER,INTENT(IN) :: m1,m2,nb
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
#endif
  END SUBROUTINE mpimat_create_blacsgrid
END MODULE m_types_mpimat
