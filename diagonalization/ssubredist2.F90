MODULE m_subredist2
CONTAINS
!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!DEC$ FREEFORM
SUBROUTINE subredist2(n,m,lda,SUB_COMM,nprow,npcol,&
     iam,ierr,nb,achi_r,asca_r,achi_c,asca_c)
#include"./cpp_double.h"
  !#include"cpp_arch.h"
  USE m_juDFT
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  ! Will compute the redistribution of the n*m-matrix
  ! Asca distributed in a block-cyclic 2-dimensional manner with
  ! block-size nb coming from ScaLAPACK
  ! to a column cyclically distributed matrix
  ! Achi used as output of chani.
  ! nprow, npcol, n,m, and nb, the blocking size for ScaLAPACK distribution
  ! will be handed to the routine.
  ! The number of processors np=nprow*npcol will be computed
  !
  ! Parameters, I/O channel numbers
  REAL,OPTIONAL     :: achi_r(*) ! Output, matrix in row cyclic-distribution
  REAL,OPTIONAL     :: asca_r(lda,*) ! Input, matrix in ScaLAPACK distribution
  COMPLEX,OPTIONAL  :: achi_c(*) ! Input, matrix in chani-distribution
  COMPLEX,OPTIONAL  :: asca_c(lda,*) ! Output, matrix in ScaLAPACK distribution

  INTEGER  :: n,m,lda   ! Global matrix sizes, local leading dimension of asca
  INTEGER  :: SUB_COMM,nprow,npcol,iam ! Communicator, 
  !                number of processor rows and columns in SUB_COMM, my rank
  INTEGER  :: nb   ! blocking size, will be computed later
  INTEGER  :: ierr  ! Error parameter to report problems
  !
  ! Local Variables
  INTEGER :: np       ! number of procs, nprow*npcol=np
  INTEGER :: myrow, mycol ! my processor coordinates
  INTEGER :: j,i,k,l, kk           ! counters
  INTEGER :: mynumcols    ! number of columns processor iam owns after redistribution
  INTEGER :: jproc, jprocrow       ! processor (row) number counter variable
  INTEGER :: hiscolssca,hisrowssca ! number of matrix columns/rows sending proc
  INTEGER :: iworldgroup           ! mpi group index for group to SUB_COMM
  INTEGER :: myrowssca,mycolssca   ! number of rows/cols processor iam has in Asca
  INTEGER, ALLOCATABLE :: irowgroup(:), icolgroup(:) ! Subgroups rows, columns
  INTEGER, ALLOCATABLE :: icommrow(:), icommcol(:) ! Communicators for Subgroups
  INTEGER, ALLOCATABLE :: iranks(:)  ! processor ranks to create subgroups
  INTEGER, ALLOCATABLE :: iwork(:,:), indx(:), isendcnts(:)  ! see comments where used
  INTEGER, ALLOCATABLE :: irecvcnts(:)  ! see comments where used
  INTEGER, ALLOCATABLE :: irecvcol(:)  ! see comments where used
  INTEGER, ALLOCATABLE :: ihavecolnum(:)  ! Which columns do I have after first alltoall
  INTEGER, ALLOCATABLE :: istartrecvbuf(:)  ! Start of cols from proc i in recvbuf
  INTEGER, ALLOCATABLE :: iinsendbuf(:)  ! Now in sendbuf for proc i
  INTEGER :: lensendbuf, lenrecvbuf ! Length of complete send- and recvbuf, first alltoall
  INTEGER :: lenrecvbuf2,lensendbuf2 ! Length of complete send- and recvbuf, second alltoall
  INTEGER :: jglob, jloc    ! Global and local column index
  INTEGER :: iglob, iloc    ! Global and local row index
  INTEGER :: hisjloc      ! Local Achi column number of global column jglob sent
  INTEGER :: numcols      !  number of columns I get in first alltoall
  INTEGER :: ierr1    ! For all the MPI calls
  !
  ! Arrays to redistribute
  !
  REAL,ALLOCATABLE :: sendbuffer_r(:)
  REAL,ALLOCATABLE :: recvbuffer_r(:)
  COMPLEX,ALLOCATABLE :: sendbuffer_c(:)
  COMPLEX,ALLOCATABLE :: recvbuffer_c(:)

  LOGICAL:: l_real
  l_real=PRESENT(achi_r)
  !
  ! Now do the computation
  !
  ierr=0
  np=nprow*npcol   ! number of processors
  mynumcols=m/np   ! my number of columns in the column-cyclic distribution
  IF (iam.LT.m-(m/np)*np) mynumcols=mynumcols+1
  myrow=iam/npcol  ! my row number in the BLACS nprow*npcol grid
  mycol=iam -(iam/npcol)*npcol  ! my column number in the BLACS nprow*npcol grid
  !
  ! Look, how many rows and columns Asca has locally
  !
  myrowssca=(n-1)/(nb*nprow)*nb+ &
       MIN(MAX(n-(n-1)/(nb*nprow)*nb*nprow-nb*myrow,0),nb)
  IF (myrowssca.GT.lda) THEN
     WRITE(6,*) 'Redist2: Wrong dimension of Asca, lda=',lda,  &
          ' myrowssca=',myrowssca
  ENDIF
  mycolssca=(m-1)/(nb*npcol)*nb+ &
       MIN(MAX(m-(m-1)/(nb*npcol)*nb*npcol-nb*mycol,0),nb)
  !
  ! Allocate arrays to build subgroups for rows and columns
  ! and create subgroups consisting of processor rows and columns
  !
  CALL MPI_COMM_GROUP(SUB_COMM,iworldgroup, Ierr)
  IF (ierr.NE.0) THEN
     WRITE(6,*) 'Mpi_comm_group Failed, ierr=',ierr
     !CALL CPP_flush(6)
     CALL juDFT_error('MPI_COMM_GROUP failed')

  ENDIF
  ALLOCATE(iranks(1:MAX(npcol,nprow)))
  ALLOCATE(irowgroup(0:nprow-1))
  ALLOCATE(icommrow(0:nprow-1))
  ALLOCATE(icolgroup(0:npcol-1))
  ALLOCATE(icommcol(0:npcol-1))
  DO k=0,nprow-1  ! Create groups of one processor row each
     DO i=1,npcol
        iranks(i)=k*npcol+i-1
     ENDDO
     CALL MPI_GROUP_INCL(iworldgroup,npcol,iranks,irowgroup(k),&
          ierr)
     IF (ierr.NE.0) THEN
        WRITE(6,*) 'MPI_GROUP_INCL(',k,') row failed, ierr=',ierr
        !CALL CPP_flush(6)
        CALL juDFT_error("MPI_GROUP_INCL failed")
     ENDIF
     CALL MPI_COMM_CREATE(SUB_COMM,irowgroup(k),&
          icommrow(k),ierr)
     IF (ierr.NE.0) THEN
        WRITE(6,*) 'MPI_COMM_CREATE(',k,') row failed, ierr=',ierr
        !CALL CPP_flush(6)
        CALL juDFT_error("MPI_COMM_CREATE failed")
     ENDIF
  ENDDO
  DO k=0,npcol-1  ! Create groups of one processor column each
     DO i=1,nprow
        iranks(i)=k+(i-1)*npcol
     ENDDO
     CALL MPI_GROUP_INCL(iworldgroup,nprow,iranks,icolgroup(k),&
          ierr)
     IF (ierr.NE.0) THEN
        WRITE(6,*) 'MPI_GROUP_INCL(',k,') col failed, ierr=',ierr
        !CALL CPP_flush(6)
        CALL juDFT_error("MPI_GROUP_INCL failed")
     ENDIF
     CALL MPI_COMM_CREATE(SUB_COMM,icolgroup(k),&
          icommcol(k),ierr)
     IF (ierr.NE.0) THEN
        WRITE(6,*) 'MPI_COMM_CREATE(',k,') col failed, ierr=',ierr
        !CALL CPP_flush(6)
        CALL juDFT_ERROR("MPI_COMM_CREATE failed")
     ENDIF
  ENDDO
  DEALLOCATE(iranks)
  !
  !  Now go through the local columns of Asca and look, to which
  !  processor column they will belong
  !
  ALLOCATE(iwork(1:mycolssca,0:MAX(npcol-1,nprow-1)))
  ! iwork(i,jproc) will tell which local ScaLAPACK column is the
  ! i-th column that will be sent to processor (myrow,modulo(jproc,npcol))
  ! in the first alltoallv
  iwork=0
  ALLOCATE(indx(0:MAX(npcol-1,nprow-1)))
  ! indx(modulo(jproc,npcol)) will indicate how many columns of processor column
  ! modulo(jproc,npcol) are already found on this processor
  indx=0  ! I have not found any column of any other processor yet
  k=0
  DO jloc=1,mycolssca,nb
     DO j=0,MIN(nb-1,mycolssca-jloc)
        jglob=k*nb*npcol+mycol*nb+J+1  ! global column jglob, local jloc+j
        jproc=MODULO(jglob-1,np) ! belonging to processor jproc
        IF (jproc == iam) THEN
           !         I will get this column on achi
           !         There are (jglob-1)/np columns before this one on achi, and
           !         It starts with element myrow*nb+1 and goes through it by steps
           !         of nb*nprow
           DO l=0,myrowssca/nb  ! Go through column jglob, store to Achi
              IF (l_real) THEN
                 DO i=1,MIN(nb,myrowssca-l*nb)
                    achi_r((jglob-1)/np*n+l*nprow*nb+myrow*nb+i)= &
                         asca_r(l*nb+i,jloc+j)
                 ENDDO
              ELSE
                 DO i=1,MIN(nb,myrowssca-l*nb)
                    achi_c((jglob-1)/np*n+l*nprow*nb+myrow*nb+i)= &
                         asca_c(l*nb+i,jloc+j)
                 ENDDO
              ENDIF
           ENDDO
        ELSE IF (MODULO(jproc,npcol) /= mycol) THEN
           !         proc modulo(jglob-1,np) is not in my processor column,
           !         but in column modulo(jproc,npcol)
           !         I will send to him in first alltoallv
           indx(MODULO(jproc,npcol))=indx(MODULO(jproc,npcol))+1
           iwork(indx(MODULO(jproc,npcol)),MODULO(jproc,npcol))= &
                jloc+j
           !
           !         I will send local column jloc+j to processor myrow,modulo(jproc,npcol)
        ENDIF
     ENDDO
     k=k+1  ! dont forget, one local block is treated
  ENDDO
  ! Now I know what to send in first alltoallv:
  ALLOCATE(isendcnts(0:MAX(npcol-1,nprow-1))) ! For alltoallv
  isendcnts=0
  lensendbuf=0 ! sendbuffer is empty
  DO k=0,npcol-1
     IF (k.NE.mycol) THEN
        isendcnts(k)=indx(k)*myrowssca
        lensendbuf=lensendbuf+isendcnts(k)
     ENDIF
  ENDDO
  IF (l_real) THEN
     ALLOCATE(sendbuffer_r(lensendbuf))
  ELSE
     ALLOCATE(sendbuffer_c(lensendbuf))
  ENDIF
  ! Put elements of asca to sendbuffer
  !
  j=0
  DO k=0,npcol-1
     IF (k.NE.mycol) THEN
        DO l=1,indx(k)
           IF (l_real) THEN
              DO i=1,myrowssca
                 sendbuffer_r(j+i)=asca_r(i,iwork(l,k))
              ENDDO
           ELSE
              DO i=1,myrowssca
                 sendbuffer_c(j+i)=asca_c(i,iwork(l,k))
              ENDDO
           END IF
           j=j+myrowssca
        ENDDO
     ENDIF
  ENDDO
  !
  ! Now look what I get in the first alltoallv
  ! Go through the other processors in my row and see what they have for me
  !
  ALLOCATE(irecvcnts(0:MAX(npcol-1,nprow-1)))
  ALLOCATE(irecvcol(1:m))
  irecvcnts=0
  irecvcol=0
  lenrecvbuf=0
  numcols=0 ! Number of columns I will get in this alltoall
  DO k=0,npcol-1 ! Now look what processors in my row send to me
     IF (k.NE.mycol) THEN  ! All except me are interesting
        hiscolssca=(m-1)/(nb*npcol)*nb+ &
             MIN(MAX(m-(m-1)/(nb*npcol)*nb*npcol-nb*k,0),nb)
        !          His number of columns in ScaLPACK distribution
        l=0
        DO jloc=1,hiscolssca,nb ! Blocks of nb columns
           DO j=0,MIN(nb-1,hiscolssca-jloc)
              jglob=l*nb*npcol+k*nb+J+1  ! global column jglob, local jloc+j
              jproc=MODULO(jglob-1,np) ! Belonging to processor jproc for achi
              IF ( MODULO(jproc,npcol).EQ.mycol) THEN
                 !               This proc is in my proc column
                 !               I will get this column in first alltoallv
                 numcols=numcols+1 ! I get one more column
                 irecvcol(numcols)=jglob ! it is global column
                 irecvcnts(k)=irecvcnts(k)+myrowssca ! more elements to rec
              ENDIF
           ENDDO
           l=l+1 ! dont forget, one local block is treated
        ENDDO
     ENDIF       ! If it is not me
     lenrecvbuf=lenrecvbuf+irecvcnts(k)
  ENDDO       ! Loop over procs in my row
  !
  ! Allocate receivebuffer
  !
  IF (l_real) THEN
     ALLOCATE(recvbuffer_r(lenrecvbuf))
  ELSE
     ALLOCATE(recvbuffer_c(lenrecvbuf))
  END IF
  ALLOCATE(iinsendbuf(0:MAX(npcol-1,nprow-1)))
  iinsendbuf(0)=0
  ALLOCATE(istartrecvbuf(0:MAX(npcol-1,nprow-1)))
  istartrecvbuf(0)=0
  DO l=0,npcol-2
     IF (l.EQ.mycol) THEN
        iinsendbuf(l+1)=iinsendbuf(l)
        istartrecvbuf(l+1)=istartrecvbuf(l)
     ELSE
        iinsendbuf(l+1)=iinsendbuf(l)+isendcnts(l)
        istartrecvbuf(l+1)=istartrecvbuf(l)+irecvcnts(l)
     ENDIF
  ENDDO
  IF(l_real) THEN
     CALL MPI_ALLTOALLV(sendbuffer_r,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_REAL,recvbuffer_r,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_REAL,icommrow(myrow),ierr)
  ELSE
     CALL MPI_ALLTOALLV(sendbuffer_c,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_COMPLEX,recvbuffer_c,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_COMPLEX,icommrow(myrow),ierr)
  ENDIF
  !
  ! After first alltoallv deallocate sendbuffer (not necessary)
  ! And set indexbuffers of sendbuffer to zero
  !
  ! deallocate(sendbuffer)
  lensendbuf2=0
  iinsendbuf(0)=0
  isendcnts=0
  !
  ! Now look how many columns of Achi the processors in my processor
  ! column have. I should have myrossca elements of all of these rows
  ! either in my Asca or in my receivebuffer
  !
  indx=m/np ! Each processor has at least m/np columns
  DO k=0,nprow-1 ! Procs in my column
     !   global processor number is k*npcol+mycol, first ones possibly one more column
     IF (k*npcol+mycol.LT.m-(m/np*np)) indx(k)=indx(k)+1
     IF (k.NE.myrow) THEN
        isendcnts(k)=indx(k)*myrowssca
        lensendbuf2=lensendbuf2+isendcnts(k)
        IF (k.LT.nprow-1) iinsendbuf(k+1)=iinsendbuf(k)+isendcnts(k)
     ELSE
        isendcnts(myrow)=0
        IF (k.LT.nprow-1) iinsendbuf(k+1)=iinsendbuf(k)
     ENDIF
  ENDDO ! End k-loop over procs in my column
  IF (lensendbuf2 > lensendbuf) THEN
     IF (l_real) THEN
        DEALLOCATE(sendbuffer_r)
        ALLOCATE(sendbuffer_r(lensendbuf2))
     ELSE
        DEALLOCATE(sendbuffer_c)
        ALLOCATE(sendbuffer_c(lensendbuf2))
     END IF

  ENDIF
  !
  ! Now look through Asca which columns will be sent now
  ! And put them to the correct place in sendbuf
  !
  k=0
  DO jloc=1,mycolssca,nb
     DO j=0,MIN(nb-1,mycolssca-jloc)
        jglob=k*nb*npcol+mycol*nb+J+1  ! global column jglob, local jloc+j
        jproc=MODULO(jglob-1,np) ! belonging to processor jproc
        !       if (jproc.eq.iam) then
        !
        !         Everything done before first alltoall, columns already on achi
        !
        !       else if (modulo(jproc,npcol).ne.mycol) then
        !
        !         proc modulo(jglob-1,np) is not in my processor column, 
        !         I sent it to him in first alltoallv
        !
        IF (MODULO(jproc,npcol) == mycol.AND.jproc /= iam) THEN
           !
           !         proc modulo(jglob-1,np) is in my column, I will send to him
           !         now in second alltoallv
           !
           jprocrow=MODULO(jglob-1,np)/npcol ! and in this row
           !         Now look which local column number for this proc column
           !         jglob will be and put it to the correct place in sendbuffer
           !         This place is iinsendbuf(k)+
           !         number of local columns of proc jprocrow,mycol before column jglob
           !         * myrowssca + i,i=1,myrowssca
           hisjloc=(jglob-1)/np+1 ! his local Achi column number
           !
           !         Store to sendbuffer
           !
           IF (l_real) THEN
              DO i=1,myrowssca
                 sendbuffer_r(iinsendbuf(jprocrow)+(hisjloc-1)*myrowssca+i) =asca_r(i,jloc+j)
              ENDDO
           ELSE
              DO i=1,myrowssca
                 sendbuffer_c(iinsendbuf(jprocrow)+(hisjloc-1)*myrowssca+i) =asca_c(i,jloc+j)
              ENDDO
           ENDIF
        ENDIF
     ENDDO
     k=k+1  ! dont forget, one local block is treated
  ENDDO
  !
  ! and through receivebuffer which columns are now there and
  ! which will be sent on
  !
  DO k=1,numcols !number of columns I got
     jglob=irecvcol(k) ! The global column number of column k in the recvbuf
     jproc=MODULO(jglob-1,np) ! The processor it will belong to is jproc
     IF (jproc.EQ.iam) THEN
        !
        !      It is my column, have to store it to achi
        !
        jloc=(jglob-1)/np+1 ! Local Achi column number
        l=0  ! Number of blocks already treated
        DO i=1,myrowssca,nb ! Go through row blocks
           IF (l_real) THEN
              DO kk=0,MIN(nb-1,myrowssca-i)
                 achi_r((jloc-1)*n+myrow*nb+l*nb*nprow+kk+1) =recvbuffer_r((k-1)*myrowssca+i+kk)
              ENDDO
           ELSE
              DO kk=0,MIN(nb-1,myrowssca-i)
                 achi_c((jloc-1)*n+myrow*nb+l*nb*nprow+kk+1) =recvbuffer_c((k-1)*myrowssca+i+kk)
              ENDDO
           ENDIF
           l=l+1 ! Next row block
        ENDDO
     ELSE  ! It is not my column, must be sent to jprocrow
        jprocrow=MODULO(jglob-1,np)/npcol ! jproc's row number '
        jloc=(jglob-1)/np+1 ! Local Achi column number for jproc
        !
        !      Store it to sendbuffer
        !
        IF (l_real) THEN
           DO i=1,myrowssca
              sendbuffer_r(iinsendbuf(jprocrow)+(jloc-1)*myrowssca+i)  &
                   =recvbuffer_r((k-1)*myrowssca+i)
           ENDDO
        ELSE
           DO i=1,myrowssca
              sendbuffer_c(iinsendbuf(jprocrow)+(jloc-1)*myrowssca+i)  &
                   =recvbuffer_c((k-1)*myrowssca+i)
           ENDDO
        ENDIF
     ENDIF
  ENDDO
  !
  !  Now look how large the next receivebuffer has to be
  !
  lenrecvbuf2=0
  irecvcnts=0
  istartrecvbuf=0
  DO k=0,nprow-1 ! Procs in my column
     IF (k.NE.myrow) THEN  ! I get from this proc
        hisrowssca= (n-1)/(nb*nprow)*nb+ &
             MIN(MAX(n-(n-1)/(nb*nprow)*nb*nprow-nb*k,0),nb) ! Proc k,mycol has
        !                                        hisrowssca rows in ScaLAPACK distribution
        !      It sends indx(myrow) columns of that length because it has elements of
        !      all of my columns
        irecvcnts(k)=indx(myrow)*hisrowssca ! I'll get irecvcnts(k) elements from '
        !                                            ( k, mycol)
        lenrecvbuf2=lenrecvbuf2+irecvcnts(k)
        IF (k.LT.nprow-1) istartrecvbuf(k+1)=  &
             istartrecvbuf(k)+irecvcnts(k)
     ELSE  ! Nothing to receive, it is me
        irecvcnts(myrow)=0
        IF (k.LT.nprow-1) istartrecvbuf(k+1)=istartrecvbuf(k)
     ENDIF
  ENDDO ! End k-loop over procs in my column
  IF (lenrecvbuf2 > lenrecvbuf) THEN
     IF (l_real) THEN
        DEALLOCATE(recvbuffer_r)
        ALLOCATE(recvbuffer_r(lenrecvbuf2))
     ELSE
        DEALLOCATE(recvbuffer_c)
        ALLOCATE(recvbuffer_c(lenrecvbuf2))
     ENDIF
  ENDIF
  !
  IF (l_real) THEN
     CALL MPI_ALLTOALLV(sendbuffer_r,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_REAL,recvbuffer_r,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_REAL,icommcol(mycol),ierr)
  ELSE
     CALL MPI_ALLTOALLV(sendbuffer_c,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_COMPLEX,recvbuffer_c,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_COMPLEX,icommcol(mycol),ierr)
  ENDIF
  !
  ! After second alltoallv deallocate sendbuffer
  ! And indexbuffers
  !
  DEALLOCATE(iwork)
  DEALLOCATE(isendcnts)
  DEALLOCATE(iinsendbuf)
  !
  ! Store the received elements to Achi
  !
  DO k=0,nprow-1 ! All except me have sent something
     IF (k.NE.myrow) THEN
        hisrowssca= (n-1)/(nb*nprow)*nb+ &
             MIN(MAX(n-(n-1)/(nb*nprow)*nb*nprow-nb*k,0),nb)
        DO j=1,indx(myrow)  ! Go through my columns, j is local Achi col num
           l=0  ! Number of blocks of this column already treated
           DO i=1,hisrowssca,nb ! Go through his rowblocks
              IF (l_real) THEN
                 DO kk=0,MIN(nb-1,hisrowssca-i)
                    achi_r((j-1)*n+k*nb+l*nb*nprow+kk+1) =recvbuffer_r(istartrecvbuf(k)+(j-1)*hisrowssca+i+kk)
                 ENDDO
              ELSE
                 DO kk=0,MIN(nb-1,hisrowssca-i)
                    achi_c((j-1)*n+k*nb+l*nb*nprow+kk+1) =recvbuffer_c(istartrecvbuf(k)+(j-1)*hisrowssca+i+kk)
                 ENDDO
              ENDIF
              l=l+1 ! Next row block
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  DEALLOCATE(indx)
  DEALLOCATE(irecvcnts)
  DEALLOCATE(istartrecvbuf)
  !
  ! Free the row and column communicators
  !
  DO k=0,nprow-1
     IF (icommrow(k).NE.MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(icommrow(k),ierr)
     ENDIF
     CALL MPI_GROUP_FREE(irowgroup(k),ierr)
  ENDDO
  DO k=0,npcol-1
     IF (icommcol(k).NE.MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(icommcol(k),ierr)
     ENDIF
     CALL MPI_GROUP_FREE(icolgroup(k),ierr)
  ENDDO
  !
  !  Deallocate the arrays to store the indices
  !
  DEALLOCATE(irecvcol)
  DEALLOCATE(irowgroup)
  DEALLOCATE(icommrow)
  DEALLOCATE(icolgroup)
  DEALLOCATE(icommcol)
  RETURN
END SUBROUTINE subredist2
END
