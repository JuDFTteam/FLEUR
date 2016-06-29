!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!DEC$ FREEFORM
subroutine subredist2(achi,n,m,asca,lda,SUB_COMM,nprow,npcol,&
                       iam,ierr,nb)
#include"./cpp_double.h"
!#include"cpp_arch.h"
use m_juDFT
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
#ifdef CPP_INVERSION
 Real     :: achi(*) ! Output, matrix in row cyclic-distribution
 Real     :: asca(lda,*) ! Input, matrix in ScaLAPACK distribution
#else
 Complex  :: achi(*) ! Input, matrix in chani-distribution
 Complex  :: asca(lda,*) ! Output, matrix in ScaLAPACK distribution
#endif
Integer  :: n,m,lda   ! Global matrix sizes, local leading dimension of asca
Integer  :: SUB_COMM,nprow,npcol,iam ! Communicator, 
!                number of processor rows and columns in SUB_COMM, my rank
Integer  :: nb   ! blocking size, will be computed later
Integer  :: ierr  ! Error parameter to report problems
!
! Local Variables
Integer :: np       ! number of procs, nprow*npcol=np
Integer :: myrow, mycol ! my processor coordinates
Integer :: j,i,k,l, kk           ! counters
Integer :: mynumcols    ! number of columns processor iam owns after redistribution
Integer :: jproc, jprocrow       ! processor (row) number counter variable
Integer :: hiscolssca,hisrowssca ! number of matrix columns/rows sending proc
Integer :: iworldgroup           ! mpi group index for group to SUB_COMM
Integer :: myrowssca,mycolssca   ! number of rows/cols processor iam has in Asca
Integer, allocatable :: irowgroup(:), icolgroup(:) ! Subgroups rows, columns
Integer, allocatable :: icommrow(:), icommcol(:) ! Communicators for Subgroups
Integer, allocatable :: iranks(:)  ! processor ranks to create subgroups
Integer, allocatable :: iwork(:,:), indx(:), isendcnts(:)  ! see comments where used
Integer, allocatable :: irecvcnts(:)  ! see comments where used
Integer, allocatable :: irecvcol(:)  ! see comments where used
Integer, allocatable :: ihavecolnum(:)  ! Which columns do I have after first alltoall
Integer, allocatable :: istartrecvbuf(:)  ! Start of cols from proc i in recvbuf
Integer, allocatable :: iinsendbuf(:)  ! Now in sendbuf for proc i
Integer :: lensendbuf, lenrecvbuf ! Length of complete send- and recvbuf, first alltoall
Integer :: lenrecvbuf2,lensendbuf2 ! Length of complete send- and recvbuf, second alltoall
Integer :: jglob, jloc    ! Global and local column index
Integer :: iglob, iloc    ! Global and local row index
Integer :: hisjloc      ! Local Achi column number of global column jglob sent
Integer :: numcols      !  number of columns I get in first alltoall
Integer :: ierr1    ! For all the MPI calls
!
! Arrays to redistribute
!
#ifdef CPP_INVERSION
 Real,allocatable :: sendbuffer(:)
 Real,allocatable :: recvbuffer(:)
#else
 Complex,allocatable :: sendbuffer(:)
 Complex,allocatable :: recvbuffer(:)
#endif
!
! Now do the computation
!
 ierr=0
 np=nprow*npcol   ! number of processors
 mynumcols=m/np   ! my number of columns in the column-cyclic distribution
 if (iam.lt.m-(m/np)*np) mynumcols=mynumcols+1
 myrow=iam/npcol  ! my row number in the BLACS nprow*npcol grid
 mycol=iam -(iam/npcol)*npcol  ! my column number in the BLACS nprow*npcol grid
!
! Look, how many rows and columns Asca has locally
!
 myrowssca=(n-1)/(nb*nprow)*nb+ &
       min(max(n-(n-1)/(nb*nprow)*nb*nprow-nb*myrow,0),nb)
 if (myrowssca.gt.lda) then
    write(6,*) 'Redist2: Wrong dimension of Asca, lda=',lda,  &
             ' myrowssca=',myrowssca
 endif
 mycolssca=(m-1)/(nb*npcol)*nb+ &
       min(max(m-(m-1)/(nb*npcol)*nb*npcol-nb*mycol,0),nb)
!
! Allocate arrays to build subgroups for rows and columns
! and create subgroups consisting of processor rows and columns
!
 Call MPI_COMM_GROUP(SUB_COMM,iworldgroup, Ierr)
 IF (ierr.ne.0) then
    write(6,*) 'Mpi_comm_group Failed, ierr=',ierr
    !CALL CPP_flush(6)
    CALL juDFT_error('MPI_COMM_GROUP failed')

 ENDIF
 allocate(iranks(1:max(npcol,nprow)))
 allocate(irowgroup(0:nprow-1))
 allocate(icommrow(0:nprow-1))
 allocate(icolgroup(0:npcol-1))
 allocate(icommcol(0:npcol-1))
 do k=0,nprow-1  ! Create groups of one processor row each
    do i=1,npcol
       iranks(i)=k*npcol+i-1
    enddo
    call MPI_GROUP_INCL(iworldgroup,npcol,iranks,irowgroup(k),&
           ierr)
    IF (ierr.ne.0) then
       write(6,*) 'MPI_GROUP_INCL(',k,') row failed, ierr=',ierr
       !CALL CPP_flush(6)
       call juDFT_error("MPI_GROUP_INCL failed")
    endif
    call MPI_COMM_CREATE(SUB_COMM,irowgroup(k),&
           icommrow(k),ierr)
    IF (ierr.ne.0) then
       write(6,*) 'MPI_COMM_CREATE(',k,') row failed, ierr=',ierr
       !CALL CPP_flush(6)
       call juDFT_error("MPI_COMM_CREATE failed")
    endif
 enddo
 do k=0,npcol-1  ! Create groups of one processor column each
    do i=1,nprow
       iranks(i)=k+(i-1)*npcol
    enddo
    call MPI_GROUP_INCL(iworldgroup,nprow,iranks,icolgroup(k),&
          ierr)
    IF (ierr.ne.0) then
       write(6,*) 'MPI_GROUP_INCL(',k,') col failed, ierr=',ierr
       !CALL CPP_flush(6)
       call juDFT_error("MPI_GROUP_INCL failed")
    endif
    call MPI_COMM_CREATE(SUB_COMM,icolgroup(k),&
         icommcol(k),ierr)
    IF (ierr.ne.0) then
       write(6,*) 'MPI_COMM_CREATE(',k,') col failed, ierr=',ierr
       !CALL CPP_flush(6)
       call juDFT_ERROR("MPI_COMM_CREATE failed")
    endif
 enddo
 deallocate(iranks)
!
!  Now go through the local columns of Asca and look, to which
!  processor column they will belong
!
 allocate(iwork(1:mycolssca,0:max(npcol-1,nprow-1)))
! iwork(i,jproc) will tell which local ScaLAPACK column is the
! i-th column that will be sent to processor (myrow,modulo(jproc,npcol))
! in the first alltoallv
 iwork=0
 allocate(indx(0:max(npcol-1,nprow-1)))
! indx(modulo(jproc,npcol)) will indicate how many columns of processor column
! modulo(jproc,npcol) are already found on this processor
 indx=0  ! I have not found any column of any other processor yet
 k=0
 Do jloc=1,mycolssca,nb
    Do j=0,min(nb-1,mycolssca-jloc)
       jglob=k*nb*npcol+mycol*nb+J+1  ! global column jglob, local jloc+j
       jproc=modulo(jglob-1,np) ! belonging to processor jproc
       if (jproc == iam) then
!         I will get this column on achi
!         There are (jglob-1)/np columns before this one on achi, and
!         It starts with element myrow*nb+1 and goes through it by steps
!         of nb*nprow
          do l=0,myrowssca/nb  ! Go through column jglob, store to Achi
             do i=1,min(nb,myrowssca-l*nb)
                achi((jglob-1)/np*n+l*nprow*nb+myrow*nb+i)= &
                   asca(l*nb+i,jloc+j)
             enddo
          enddo
       else if (modulo(jproc,npcol) /= mycol) then
!         proc modulo(jglob-1,np) is not in my processor column,
!         but in column modulo(jproc,npcol)
!         I will send to him in first alltoallv
          indx(modulo(jproc,npcol))=indx(modulo(jproc,npcol))+1
          iwork(indx(modulo(jproc,npcol)),modulo(jproc,npcol))= &
             jloc+j
!
!         I will send local column jloc+j to processor myrow,modulo(jproc,npcol)
       endif
    Enddo
    k=k+1  ! dont forget, one local block is treated
 Enddo
! Now I know what to send in first alltoallv:
 allocate(isendcnts(0:max(npcol-1,nprow-1))) ! For alltoallv
 isendcnts=0
 lensendbuf=0 ! sendbuffer is empty
 Do k=0,npcol-1
    If (k.ne.mycol) then
       isendcnts(k)=indx(k)*myrowssca
       lensendbuf=lensendbuf+isendcnts(k)
    Endif
 Enddo
 allocate(sendbuffer(lensendbuf))
!
! Put elements of asca to sendbuffer
!
 j=0
 Do k=0,npcol-1
    If (k.ne.mycol) then
       Do l=1,indx(k)
          Do i=1,myrowssca
             sendbuffer(j+i)=asca(i,iwork(l,k))
          Enddo
          j=j+myrowssca
       Enddo
    Endif
 Enddo
!
! Now look what I get in the first alltoallv
! Go through the other processors in my row and see what they have for me
!
 allocate(irecvcnts(0:max(npcol-1,nprow-1)))
 allocate(irecvcol(1:m))
 irecvcnts=0
 irecvcol=0
 lenrecvbuf=0
 numcols=0 ! Number of columns I will get in this alltoall
 do k=0,npcol-1 ! Now look what processors in my row send to me
    if (k.ne.mycol) then  ! All except me are interesting
       hiscolssca=(m-1)/(nb*npcol)*nb+ &
           min(max(m-(m-1)/(nb*npcol)*nb*npcol-nb*k,0),nb)
!          His number of columns in ScaLPACK distribution
       l=0
       Do jloc=1,hiscolssca,nb ! Blocks of nb columns
          Do j=0,min(nb-1,hiscolssca-jloc)
             jglob=l*nb*npcol+k*nb+J+1  ! global column jglob, local jloc+j
             jproc=modulo(jglob-1,np) ! Belonging to processor jproc for achi
             IF ( modulo(jproc,npcol).eq.mycol) then
!               This proc is in my proc column
!               I will get this column in first alltoallv
                numcols=numcols+1 ! I get one more column
                irecvcol(numcols)=jglob ! it is global column
                irecvcnts(k)=irecvcnts(k)+myrowssca ! more elements to rec
             Endif
          Enddo
          l=l+1 ! dont forget, one local block is treated
       Enddo
    endif       ! If it is not me
    lenrecvbuf=lenrecvbuf+irecvcnts(k)
 enddo       ! Loop over procs in my row
!
! Allocate receivebuffer
!
 allocate(recvbuffer(lenrecvbuf))
 allocate(iinsendbuf(0:max(npcol-1,nprow-1)))
 iinsendbuf(0)=0
 allocate(istartrecvbuf(0:max(npcol-1,nprow-1)))
 istartrecvbuf(0)=0
 do l=0,npcol-2
    if (l.eq.mycol) then
       iinsendbuf(l+1)=iinsendbuf(l)
       istartrecvbuf(l+1)=istartrecvbuf(l)
    else
       iinsendbuf(l+1)=iinsendbuf(l)+isendcnts(l)
       istartrecvbuf(l+1)=istartrecvbuf(l)+irecvcnts(l)
    endif
 enddo
#ifdef CPP_INVERSION
 Call MPI_ALLTOALLV(sendbuffer,isendcnts,iinsendbuf,&
       CPP_MPI_TYP_REAL,recvbuffer,irecvcnts,istartrecvbuf,&
       CPP_MPI_TYP_REAL,icommrow(myrow),ierr)
#else
 Call MPI_ALLTOALLV(sendbuffer,isendcnts,iinsendbuf,&
       CPP_MPI_TYP_COMPLEX,recvbuffer,irecvcnts,istartrecvbuf,&
       CPP_MPI_TYP_COMPLEX,icommrow(myrow),ierr)
#endif
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
 do k=0,nprow-1 ! Procs in my column
!   global processor number is k*npcol+mycol, first ones possibly one more column
    if (k*npcol+mycol.lt.m-(m/np*np)) indx(k)=indx(k)+1
    if (k.ne.myrow) then
       isendcnts(k)=indx(k)*myrowssca
       lensendbuf2=lensendbuf2+isendcnts(k)
       if (k.lt.nprow-1) iinsendbuf(k+1)=iinsendbuf(k)+isendcnts(k)
    else
       isendcnts(myrow)=0
       if (k.lt.nprow-1) iinsendbuf(k+1)=iinsendbuf(k)
    endif
 enddo ! End k-loop over procs in my column
 If (lensendbuf2 > lensendbuf) then
    deallocate(sendbuffer)
    allocate(sendbuffer(lensendbuf2))
 Endif
!
! Now look through Asca which columns will be sent now
! And put them to the correct place in sendbuf
!
 k=0
 Do jloc=1,mycolssca,nb
    Do j=0,min(nb-1,mycolssca-jloc)
       jglob=k*nb*npcol+mycol*nb+J+1  ! global column jglob, local jloc+j
       jproc=modulo(jglob-1,np) ! belonging to processor jproc
!       if (jproc.eq.iam) then
!
!         Everything done before first alltoall, columns already on achi
!
!       else if (modulo(jproc,npcol).ne.mycol) then
!
!         proc modulo(jglob-1,np) is not in my processor column, 
!         I sent it to him in first alltoallv
!
       If (modulo(jproc,npcol) == mycol.and.jproc /= iam) then
!
!         proc modulo(jglob-1,np) is in my column, I will send to him
!         now in second alltoallv
!
          jprocrow=modulo(jglob-1,np)/npcol ! and in this row
!         Now look which local column number for this proc column
!         jglob will be and put it to the correct place in sendbuffer
!         This place is iinsendbuf(k)+
!         number of local columns of proc jprocrow,mycol before column jglob
!         * myrowssca + i,i=1,myrowssca
          hisjloc=(jglob-1)/np+1 ! his local Achi column number
!
!         Store to sendbuffer
!
          do i=1,myrowssca
             sendbuffer(iinsendbuf(jprocrow)+(hisjloc-1)*myrowssca+i) &
               =asca(i,jloc+j)
          enddo
       endif
    Enddo
    k=k+1  ! dont forget, one local block is treated
 Enddo
!
! and through receivebuffer which columns are now there and
! which will be sent on
!
 Do k=1,numcols !number of columns I got
    jglob=irecvcol(k) ! The global column number of column k in the recvbuf
    jproc=modulo(jglob-1,np) ! The processor it will belong to is jproc
    if (jproc.eq.iam) then
!
!      It is my column, have to store it to achi
!
       jloc=(jglob-1)/np+1 ! Local Achi column number
       l=0  ! Number of blocks already treated
       Do i=1,myrowssca,nb ! Go through row blocks
          Do kk=0,min(nb-1,myrowssca-i)
             achi((jloc-1)*n+myrow*nb+l*nb*nprow+kk+1)  &
                 =recvbuffer((k-1)*myrowssca+i+kk)
          Enddo
          l=l+1 ! Next row block
       Enddo
    else  ! It is not my column, must be sent to jprocrow
       jprocrow=modulo(jglob-1,np)/npcol ! jproc's row number '
       jloc=(jglob-1)/np+1 ! Local Achi column number for jproc
!
!      Store it to sendbuffer
!
       Do i=1,myrowssca
          sendbuffer(iinsendbuf(jprocrow)+(jloc-1)*myrowssca+i)  &
                 =recvbuffer((k-1)*myrowssca+i)
       Enddo
    endif
 Enddo
!
!  Now look how large the next receivebuffer has to be
!
 lenrecvbuf2=0
 irecvcnts=0
 istartrecvbuf=0
 do k=0,nprow-1 ! Procs in my column
    if (k.ne.myrow) then  ! I get from this proc
       hisrowssca= (n-1)/(nb*nprow)*nb+ &
          min(max(n-(n-1)/(nb*nprow)*nb*nprow-nb*k,0),nb) ! Proc k,mycol has
!                                        hisrowssca rows in ScaLAPACK distribution
!      It sends indx(myrow) columns of that length because it has elements of
!      all of my columns
       irecvcnts(k)=indx(myrow)*hisrowssca ! I'll get irecvcnts(k) elements from '
!                                            ( k, mycol)
       lenrecvbuf2=lenrecvbuf2+irecvcnts(k)
       if (k.lt.nprow-1) istartrecvbuf(k+1)=  &
           istartrecvbuf(k)+irecvcnts(k)
    else  ! Nothing to receive, it is me
       irecvcnts(myrow)=0
       if (k.lt.nprow-1) istartrecvbuf(k+1)=istartrecvbuf(k)
    endif
 enddo ! End k-loop over procs in my column
 If (lenrecvbuf2 > lenrecvbuf) then
    deallocate(recvbuffer)
    allocate(recvbuffer(lenrecvbuf2))
 Endif
!
#ifdef CPP_INVERSION
 Call MPI_ALLTOALLV(sendbuffer,isendcnts,iinsendbuf,&
       CPP_MPI_TYP_REAL,recvbuffer,irecvcnts,istartrecvbuf,&
       CPP_MPI_TYP_REAL,icommcol(mycol),ierr)
#else
 Call MPI_ALLTOALLV(sendbuffer,isendcnts,iinsendbuf,&
       CPP_MPI_TYP_COMPLEX,recvbuffer,irecvcnts,istartrecvbuf,&
       CPP_MPI_TYP_COMPLEX,icommcol(mycol),ierr)
#endif
!
! After second alltoallv deallocate sendbuffer
! And indexbuffers
!
 deallocate(sendbuffer)
 deallocate(iwork)
 deallocate(isendcnts)
 deallocate(iinsendbuf)
!
! Store the received elements to Achi
!
 do k=0,nprow-1 ! All except me have sent something
    If (k.ne.myrow) then
       hisrowssca= (n-1)/(nb*nprow)*nb+ &
          min(max(n-(n-1)/(nb*nprow)*nb*nprow-nb*k,0),nb)
       Do j=1,indx(myrow)  ! Go through my columns, j is local Achi col num
          l=0  ! Number of blocks of this column already treated
          Do i=1,hisrowssca,nb ! Go through his rowblocks
             Do kk=0,min(nb-1,hisrowssca-i)
                achi((j-1)*n+k*nb+l*nb*nprow+kk+1)  &
                  =recvbuffer(istartrecvbuf(k)+(j-1)*hisrowssca+i+kk)
             Enddo
             l=l+1 ! Next row block
          Enddo
       Enddo
    Endif
 Enddo
 deallocate(indx)
 deallocate(irecvcnts)
 deallocate(istartrecvbuf)
!
! Free the row and column communicators
!
 do k=0,nprow-1
    if (icommrow(k).ne.MPI_COMM_NULL) then
       call MPI_COMM_FREE(icommrow(k),ierr)
    endif
    call MPI_GROUP_FREE(irowgroup(k),ierr)
 enddo
 do k=0,npcol-1
    if (icommcol(k).ne.MPI_COMM_NULL) then
       call MPI_COMM_FREE(icommcol(k),ierr)
    endif
    call MPI_GROUP_FREE(icolgroup(k),ierr)
 enddo
!
!  Deallocate the arrays to store the indices
!
 deallocate(recvbuffer)
 deallocate(irecvcol)
 deallocate(irowgroup)
 deallocate(icommrow)
 deallocate(icolgroup)
 deallocate(icommcol)
return
end subroutine subredist2
