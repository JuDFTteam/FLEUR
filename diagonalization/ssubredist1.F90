!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!DEC$ FREEFORM
SUBROUTINE subredist1(achi,n,asca,lda,SUB_COMM,nprow,npcol,iam,ierr,nb)
use m_juDFT
#include"./cpp_double.h"
!#include"cpp_arch.h"

IMPLICIT NONE
INCLUDE 'mpif.h'

! Will compute the redistribution of the half-colums processor iam owns
! in the beginning when an n*n-matrix is distributed in the way
! used in the beginning of chani:
! Processor iam owns the upper half including the diagonal of columns
! iam+(k*np)+1, k=0,...n/np[+1]
! stored contigously in a one-dimensional array achi
! to the matrix asca which is distributed to a
! logical two-dimensional ScaLAPACK-grid with blocking size nb
! For testing reasons nprow, npcol, n, and nb, the blocking size 
! for ScaLAPACK distribution
! will be handed to the routine. Later on, nb will be computed from 
! nprow, npcol, and n.
! The number of processors np=nprow*npcol will be computed
!
! Parameters, I/O channel numbers
#ifdef CPP_INVERSION
 Real     :: achi(*) ! Input, matrix in chani-distribution
 Real     :: asca(lda,*) ! Output, matrix in ScaLAPACK distribution
#else
 Complex  :: achi(*) ! Input, matrix in chani-distribution
 Complex  :: asca(lda,*) ! Output, matrix in ScaLAPACK distribution
#endif
Integer  :: n,lda   ! Global matrix size, local leading dimension of asca
Integer  :: SUB_COMM,nprow,npcol,iam ! Communicator, 
!                number of processor rows and columns in SUB_COMM, my rank
Integer  :: nb   ! blocking size, will be computed later
Integer  :: ierr  ! Error parameter to report problems
!
! Local Variables
Integer :: np       ! number of procs, nprow*npcol=np
Integer :: me, myrow, mycol ! my processor coordinates
Integer :: j,i,k,l, kk           ! counters
Integer :: mynumcols             ! number of columns processor iam owns in the beginning
Integer :: jproccol, jprocrow    ! processor column/row number counter variable
Integer :: isndproc, him         ! processor in my column sending to me
Integer :: hisnumcols            ! number of matrix columns proc isndproc owns
Integer :: iworldgroup           ! mpi group index for group to SUB_COMM
Integer :: myrowssca,mycolssca   ! number of rows/cols processor iam gets in Asca
Integer, allocatable :: irowgroup(:), icolgroup(:) ! Subgroups rows, columns
Integer, allocatable :: icommrow(:), icommcol(:) ! Communicators for Subgroups
Integer, allocatable :: iranks(:)  ! processor ranks to create subgroups
Integer, allocatable :: glcolnum(:,:)  ! Global column number of the columns I get
Integer :: jstart, jend, npfac   ! see comments where used
Integer, allocatable :: iwork(:,:), indx(:), isendcnts(:)  ! see comments where used
Integer, allocatable :: irecvcnts(:)  ! see comments where used
Integer, allocatable :: irecvcol(:)  ! see comments where used
!Integer, allocatable :: ihavecolnum(:)  ! Which columns do I have after first alltoall
Integer, allocatable :: istartrecvbuf(:)  ! Start of cols from proc i in recvbuf
Integer, allocatable :: iinsendbuf(:)  ! Now in sendbuf for proc i
Integer :: nsq      ! Integer sqrt of np, to compute grid
Integer :: lensendbuf, lenrecvbuf ! Length of complete send- and recvbuf, first send
Integer :: lensendbuf2, lenrecvbuf2 ! Length of complete send- and recvbuf, second send
Integer :: jglob, jloc    ! Global and local column index
Integer :: lglob, lloc    ! Global and local row index
!Integer :: icntcol    ! Number of columns I own after first alltoall
Integer :: indxrecvbuf    ! Where I am in receivebuffer of first alltoall
Integer :: ierr1    ! For all the MPI calls
Integer :: subsize,mysubrank    ! Test of comm_create
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
 np=nprow*npcol   ! number of processors
 mynumcols=n/np   ! my number of columns in the chani-distribution
 if (iam < n-(n/np)*np) mynumcols=mynumcols+1
 myrow=iam/npcol  ! my row number in the BLACS nprow*npcol grid
 mycol=iam -(iam/npcol)*npcol  ! my column number in the BLACS nprow*npcol grid
 me = iam+1       ! makes life easyer
!
!  Compute how many elements of Achi this processor owns (Never needed)
!
! j=0
! Do i=me,n,np ! Go through my columns of Achi
!    j=j+i ! Number of elements I have in Achi
! Enddo
!
! Get Group handle for whole group
!
 Call MPI_COMM_GROUP(SUB_COMM,iworldgroup, Ierr)
 IF (ierr /= 0) then
    write(6,*) 'MPI_COMM_GROUP failed, ierr=',ierr
!    CALL CPP_flush(6)
    call juDFT_error("MPI_COMM_GROUP failed")
 ENDIF
!
! Allocate arrays to build subgroups for rows and columns
! and create subgroups consisting of processor rows and columns
!
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
    IF (ierr /= 0) then
       write(6,*) 'MPI_GROUP_INCL(',k,') row failed, ierr=',ierr
!       CALL CPP_flush(6)
       call juDFT_error("MPI_COMM_INCL failed")
    endif
    call MPI_COMM_CREATE(SUB_COMM,irowgroup(k),&
           icommrow(k),ierr)
    IF (ierr /= 0) then
       write(6,*) 'MPI_COMM_CREATE(',k,') row failed, ierr=',ierr
!       CALL CPP_flush(6)
       call juDFT_error("MPI_COMM_CREATE failed")
    endif
 enddo
 do k=0,npcol-1  ! Create groups of one processor column each
    do i=1,nprow
       iranks(i)=k+(i-1)*npcol
    enddo
    call MPI_GROUP_INCL(iworldgroup,nprow,iranks,icolgroup(k),&
          ierr)
    IF (ierr /= 0) then
       write(6,*) 'MPI_GROUP_INCL(',k,') col failed, ierr=',ierr
       !CALL CPP_flush(6)
       call juDFT_error("MPI_COMM_INCL failed")
    endif
    call MPI_COMM_CREATE(SUB_COMM,icolgroup(k),&
         icommcol(k),ierr)
    IF (ierr /= 0) then
       write(6,*) 'MPI_COMM_CREATE(',k,') col failed, ierr=',ierr
       !CALL CPP_flush(6)
       call juDFT_error("MPI_COMM_CREATE failed")
    endif
 enddo
   deallocate(iranks)
!
!
! Allocate the arrays to store the indices
!
 allocate(indx(0:npcol-1))
 allocate(isendcnts(0:max(npcol-1,nprow-1)))
 allocate(iinsendbuf(0:max(npcol-1,nprow-1)))
! indx(jproccol) will indicate how many columns of processor column
! jproccol are already found on this processor
 indx=0  ! I have not found any column of any other processor yet
 isendcnts=0  ! The length is 0 thus
 allocate(iwork(1:2*mynumcols,0:npcol-1))
 iwork=0
!  allocate(ihavecolnum( ((n/nb+1)/npcol+1)*nb ))
!  ihavecolnum=0
! iwork(2l-1,jproccol) will show the starting and iwork(2l,jproccol)
! the end-indices for the columns of
! processor column jproccol in the one-dimensional array of processor iam
! Example: processor iam=me-1 has the upper parts of columns me+k*np, k=0,..,n/np[+1]
! contigously in its local one-dimensional array. Thus elements
! 1,..,me belong to global column jglob=me, elements me+1,...,2*me+np belong to
! global column jglob=me+np ...
! Global column jglob belongs to processor column 
! jproccol= modulo(jglob-1,nb*npcol)/nb
! Thus iwork(1,modulo(me-1,nb*npcol)/nb)=1
! iwork(2,modulo(me-1,nb*npcol)/nb)=me
! indx(modulo(me-1,nb*npcol)/nb)=indx(modulo(me-1,nb*npcol)/nb)+1
! iwork(2*indx(modulo(me+np-1,nb*npcol)/nb)+1,modulo(me+np-1,nb*npcol)/nb)=me+1
! iwork(2*indx(modulo(me+np-1,nb*npcol)/nb)+2,modulo(me+np-1,nb*npcol)/nb)=2*me+np
!
 jend=0  ! jstart will be jend+1
 npfac=0 ! my first elements belong to a column shorter or equal np
 do j=1,mynumcols  ! Find out who my columns belong to
    jstart=jend+1  ! Start point in achi 
    jend=jend+me+npfac*np ! End point in achi
    jglob=me+npfac*np  ! global column number
    npfac=npfac+1  ! next column is np elements longer
    jproccol=modulo(jglob-1,npcol*nb)/nb ! Belongs to proc column
    iwork(2*indx(jproccol)+1,jproccol)=jstart
    iwork(2*indx(jproccol)+2,jproccol)=jend
    isendcnts(jproccol)=isendcnts(jproccol)+jend-jstart+1
!   Send isendcnts(jproccol) elements to (myrow,jproccol)
    indx(jproccol)=indx(jproccol)+1 ! Send indx(jproccol) columns to him
 enddo  ! End j-loop over my columns
! Count the number of columns and elements I have to send to which processor column
! Here I only want to compute lensendbuf, thus I don't need to store
! anything to sendbuffer because sendbuffer is not yet allocated!!!
 lensendbuf=0 ! Length of sendbuffer
 do j=0,npcol-1  
    if (j /= mycol) lensendbuf=lensendbuf+isendcnts(j)
 enddo
 allocate(sendbuffer(lensendbuf))
!
!  Now I have to put the elements to be sent to sendbuffer !!
!
 lensendbuf=0
 iinsendbuf(0)=0 ! Starting addresses in sendbuffer (See MPI_ALLTOALLV)
 do j=0,npcol-1  
    if (indx(j)==0) then
       if (isendcnts(j) /= 0) write(6,*) iam, &
               ' Fehler, isendcnts nicht 0',isendcnts(j)
       if (j < npcol-1) iinsendbuf(j+1)=iinsendbuf(j)
    else
       if (j /= mycol) then
          do i=1,indx(j)
             do kk=iwork(2*i-1,j),iwork(2*i,j)
                 sendbuffer(lensendbuf+1+kk-iwork(2*i-1,j))=&
                   Achi(kk)
             enddo
             lensendbuf=lensendbuf+iwork(2*i,j)-iwork(2*i-1,j)+1
          enddo
       endif
       if (j == mycol) isendcnts(j)=0 ! must be 0 for mpi_alltoallv
       if (j < npcol-1)  iinsendbuf(j+1)=iinsendbuf(j)+isendcnts(j)
    endif
 enddo
 allocate(irecvcnts(0:max(npcol-1,nprow-1))) ! Number of elements I will get
 allocate(irecvcol(0:npcol-1)) ! Number of columns I will get
 allocate(glcolnum(n/np+1,0:npcol-1)) ! Global column numbers
 allocate(istartrecvbuf(0:max(npcol-1,nprow-1))) ! Start in receivebuffer
 irecvcnts=0
 irecvcol=0
 glcolnum=0
 lenrecvbuf=0
 istartrecvbuf(0)=0
 do k=0,npcol-1 ! Now look what processors in my row send to me
    isndproc=myrow*npcol+k ! Number of processor (myrow,k)
    if (k /= mycol) then  ! All except me are interesting
       him=isndproc+1  ! make life easier
       hisnumcols=n/np   ! his number of columns in the chani-distribution
       if (isndproc < n-(n/np)*np) hisnumcols=hisnumcols+1
       jend=0  ! jstart will be jend+1
       npfac=0 ! his first elements belong to a column shorter or equal np
       do j=1,hisnumcols  ! Look through his columns
          jstart=jend+1 ! Starting in his part of Achi
          jend=jend+him+npfac*np ! Ending in his part of Achi
          jglob=him+npfac*np ! Global column number
          npfac=npfac+1  ! next column is np elements longer
          jproccol=modulo(jglob-1,npcol*nb)/nb ! Will belong to proc column
          if (jproccol==mycol) then ! I will get this column
             irecvcol(k)=irecvcol(k)+1 ! I get one more column from (myrow,k)
             irecvcnts(k)=irecvcnts(k)+jend-jstart+1
             glcolnum(irecvcol(k),k)=jglob ! It is global column number
          endif  ! I get this column
       enddo       ! j=1,hisnumcols
       lenrecvbuf=lenrecvbuf+irecvcnts(k)
       if (k < npcol-1) istartrecvbuf(k+1) = &
                 istartrecvbuf(k)+irecvcnts(k)  
!         what I got from k starts at istartrecvbuf(k)+1 in the recvbuf
    else ! it is me, start in recvbuf for the next is same as for this 
       if (k < npcol-1) istartrecvbuf(k+1)= istartrecvbuf(k)
    endif       ! If it is not me
 enddo       ! Loop over procs in my row
 allocate(recvbuffer(lenrecvbuf))
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
!  Now look to Asca to put the elements of Achi or receivebuffer to
!  See later, where the number of elements of column jglob this
!  processor will get is computed. Jglob there plays the same role as
!  n here.
!
 myrowssca=(n-1)/(nb*nprow)*nb+ &
        min(max(n-(n-1)/(nb*nprow)*nb*nprow-nb*myrow,0),nb)
!  My number of rows in ScaLAPACK distribution
 if (myrowssca > lda) &
    write(6,*) 'Redist1: Wrong dimension of Asca, lda=',lda,  &
             ' myrowssca=',myrowssca
! My number of columns in ScaLAPACK distribution
 mycolssca=(n-1)/(nb*npcol)*nb+ &
        min(max(n-(n-1)/(nb*npcol)*nb*npcol-nb*mycol,0),nb)
!
! See computing jloc from jglob
!
 istartrecvbuf=istartrecvbuf+1 ! Makes live easier
!  
! Now look what this processor should have after the first alltoall
! along its processor row
!   
! First those columns it already had, they are still on Achi
!
 isendcnts=0 ! Number of elements alltoall to send along column
! icntcol=0
 lensendbuf2=0  ! New sendbuffer
 Do k=1,indx(mycol)    ! I had indx(mycol) columns for my processor column
!   I have to store it to Asca or sent in next Alltoallv
    jglob = iwork(2*k,mycol)-iwork(2*k-1,mycol)+1 !global column index
    jloc=jglob-(jglob-1)/(nb*npcol)*nb*(npcol-1)-nb*mycol 
!   jloc is local column index for ScaLAPACK distribution
!   Achtung, hier stand urspruenglich 
!   jloc= jglob - jglob/(nb*npcol)*nb*(npcol-1)-nb*mycol, aber das
!   fuehrte im Fall der letzten Spalte eines "sweeps" zu einem Fehler
!   Gedanke, der diese Formel ergibt:
!   Man treagt von der globalen Spaltennummer jeweils einen kompletten
!   "sweep" von nb*npcol Spalten ab, also hat man
!   jglob/(nb*npcol) mal jedem Prozessor dieser Zeile einen Spalten
!   block gegeben, also habe auch ich lokal jglob/(nb*npcol) ganze
!   Spaltenbloecke, es beleiben
!   jglob - jglob/(nb*npcol)*nb*npcol Spalten uebrig, die noch
!   zu verteilen sind, davon haben die mycol Prozessoren, die vor
!   mir in der Reihe sind jeweils einen ganzen Block,
!   bleiben also jglob - jglob/(nb*npcol)*nb*npcol - nb*mycol
!   Elemente uebrig.
!   Damit bin ich lokal bei folgender Spalte:
!   jglob/(nb*npcol) ganze "sweeps", also lokale Bloecke bis zu dieser
!   Spalte, das sind also jglob/(nb*npcol)*nb Spalten, dann noch die Restspalten,
!   Also jglob - jglob/(nb*npcol)*nb*npcol - nb*mycol Spalten.
!   Das geht gut, ausser bei der letzten Spalte eines "sweeps".
!   Hier ergibt sich jglob - jglob/(nb*npcol)*nb*npcol=0, also keine Spalten
!   uebrig. Diese letzte Spalte hat aber Prozessor (myrow, npcol-1), also
!   mycol=npcol-1, das ist meistens groesser als 0, also werden
!   nb*(npcol-1) Spalten vom Rest abgezogen, dieser Rest ist aber =0, also
!   entsteht eine negative Zahl, die zu den Spalten von den ganzen Bloecken
!   addiert wird. Eigentlich darf aber nichts mehr addiert oder subtrahiert
!   Werden, ich habe genau jglob/(nb*npcol)*nb Spalten.
!   Daher muss die Zahl der kompletten "sweeps" berechnet werden, die vor
!   dieser Spalte jglob schon fertig verteilt sind, dies ist dieselbe Anzahl
!   ausser in dem einen Fall, in dem es gerade aufgeht.
!   Wenn man dann jglob - (jglob-1)/(nb*npcol)*nb*npcol berechnet, kommt
!   genau nb*npcol, also ein ganzer "sweep" als Rest heraus. 
!   Zieht man jetzt die Spalten ab, die die Prozessoren vor mir,
!   also vor (myrow, npcol-1) haben, so erhaelt man 
!   nb*npcol-nb*mycol=nb*(npcol-mycol)=nb*(npcol-(npcol-1))=nb als lokalen Index
!   fuer den Rest, und ich bekomme als lokalen Spaltenindex
!   (Anzahl kompletter "sweeps" bis jglob-1)*nb+nb
!   =(jglob-1)/(nb*npcol)*nb+nb, und da in diesem speziellen Fall, und nur
!   in diesem Fall gilt (jglob-1)/(nb*npcol)=jglob/(nb*npcol)-1
!   erhalte ich jloc=jglob/(nb*npcol)*nb, und das stimmt, denn es ist die
!   letzte Spalte im Block jglob/(nb*npcol)
    lloc=1 ! local row index for ScaLAPACK distribution
!    icntcol= icntcol+1
!    ihavecolnum(icntcol)=jglob
    do j=1,jglob,nb
!      Go through column jglob by row blocks
!      The row block treated now belongs to processor row
!      mod(j/nb,nprow)
       if (mod(j/nb,nprow)==myrow) then ! Belonging to my processor row 
!                                         = I get it, store it to Asca
          do lglob=j,min(jglob,j+nb-1) ! This row block, no more than jglob
             Asca(lloc,jloc)=Achi(iwork(2*k-1,mycol)+lglob-1)
             lloc=lloc+1
          enddo
       else ! send to processor row mod(j/nb,nprow)
          isendcnts(mod(j/nb,nprow))=isendcnts(mod(j/nb,nprow))&
               + min(nb,jglob-j+1) ! Compute only number of elements
       endif
    enddo
 enddo
!
!  Now processing the received columns, processor by processor
!
! icntcol=indx(mycol)
 indxrecvbuf=0  ! Go through receivebuffer element by element
 do l=0,npcol-1
    do k=1,irecvcol(l) !  If l=mycol irecvcol=0, nothing will be done
       jglob=glcolnum(k,l) ! Received column is global column number jglob
       jloc=jglob-(jglob-1)/(nb*npcol)*nb*(npcol-1)-nb*mycol
!      It will be local column number jloc of Asca
!       icntcol = icntcol+1 ! Number of columns I have now
!       ihavecolnum(icntcol)=jglob ! Global indices of columns I have
       lloc=1 ! Now go through column by row blocks
       do j=1,jglob,nb
!         Go through column jglob by row blocks
!         The row block treated now belongs to processor row
!         mod(j/nb,nprow)
          if (mod(j/nb,nprow)==myrow) then  ! my row block, store it to Asca
             do lglob=j,min(jglob,j+nb-1) ! This row block, no more than jglob
                indxrecvbuf=indxrecvbuf+1
                Asca(lloc,jloc)=recvbuffer(indxrecvbuf)
                lloc=lloc+1
             enddo
          else  ! not my row block send to processor row mod(j/nb,nprow)
!                 min(nb,jglob-j+1) more elements
             isendcnts(mod(j/nb,nprow))=   &
                  isendcnts(mod(j/nb,nprow))+ min(nb,jglob-j+1)
             indxrecvbuf=indxrecvbuf+min(nb,jglob-j+1)
          endif ! who's row block
       enddo ! End going through column jglob
    enddo ! End going through columns received from processor myrow,k
 enddo ! End going through all processors in my row
! Deallocate some workspace no longer needed
 deallocate(glcolnum)
 deallocate(irecvcol)
!  deallocate(ihavecolnum)
!
 iinsendbuf(0)=1
 Do k=0,nprow-1
    If (k /= myrow) then ! I will send isendcnts(k )elements to processor
!                          k,mycol
       lensendbuf2=lensendbuf2+isendcnts(k)
       if (k < nprow-1) iinsendbuf(k+1)=iinsendbuf(k)+isendcnts(k)
    Else  ! It is my row, don't have to send
       if (k < nprow-1) iinsendbuf(k+1)=iinsendbuf(k)
    Endif
 Enddo
 if (lensendbuf2 > lensendbuf) then
    deallocate(sendbuffer)
    Allocate(sendbuffer(lensendbuf2))
 endif
!
!  Now look through all columns I need and see whether I have them or get them
!
 irecvcnts=0 ! Count for next receive
 lenrecvbuf2=0  ! New receivebuffer 
 jstart=0    ! No column from Achi stored to sendbuf
 do k=mycol*nb+1,n,npcol*nb ! k is the first element in a column block
    do jglob=k,min(k+nb-1,n) ! global column index
!      the processor originally having this column is 
!      processor modulo(jglob-1,np)
       if (modulo(jglob-1,np)==iam) then ! It was my column in the beginning
          jstart=jstart+1
          do j=1,jglob,nb ! Row block starting with j
!            Go through column jglob by row blocks
!            The row block treated now belongs to processor row
!            mod(j/nb,nprow)
             if (mod(j/nb,nprow) /= myrow) then  ! not my row block
!               Now find out where it is in Achi
                Do l=j,min(j+nb-1,jglob)
                   sendbuffer(iinsendbuf(mod(j/nb,nprow)))=&
                       Achi(iwork(2*jstart-1,mycol)+l-1)
                   iinsendbuf(mod(j/nb,nprow))= &
                      iinsendbuf(mod(j/nb,nprow))+1
                enddo
             endif ! who's row block
          enddo ! End going through column jglob
       else if (modulo(jglob-1,np)/npcol==myrow) then
          jproccol=modulo(modulo(jglob-1,np),npcol)
!         modulo(jglob-1,np) is the number of the processor who had
!         this column in the beginning
!         This processor is in my processor row, if its number divided by
!         npcol is the same as my number divided by npcol, and this is
!         myrow
!         processor number divided by npcol = processor row number
!         processor number modulo npcol = processor column number
!         It is in my receive buffer
          do j=1,jglob,nb
!            Go through column jglob by row blocks
!            The row block treated now belongs to processor row
!            mod(j/nb,nprow)
             if (mod(j/nb,nprow)==myrow) then  ! my row block
!               stored the elements from recvbuf to Asca
                istartrecvbuf(jproccol)=istartrecvbuf(jproccol)+&
                  min(nb,jglob-j+1)
             else  ! not my row block, store to sendbuffer for row mod(j/nb,nprow)
!               Now find out where it is in recvbuffer
                Do l=j,min(j+nb-1,jglob)
                   sendbuffer(iinsendbuf(mod(j/nb,nprow)))=&
                       recvbuffer(istartrecvbuf(jproccol))
                   iinsendbuf(mod(j/nb,nprow))= &
                       iinsendbuf(mod(j/nb,nprow))+1
                   istartrecvbuf(jproccol)= &
                       istartrecvbuf(jproccol)+1
                enddo
             endif ! who's row block
          enddo ! End going through column jglob
       else  ! I get it in this alltoall
          jprocrow=modulo(jglob-1,np)/npcol ! processor row where it is
!         The number of elements I have to receive for this column is:
!         The number of complete "sweeps" of nb*nprow elements of this
!         column are (jglob-1)/(nb*nprow), that is, each processor in my
!         processor column gets at least (jglob-1)/(nb*nprow)*nb elements
!         of column jglob.
!         The rest is jglob-(jglob-1)/(nb*nprow)*nb*nprow elements.
!         Out of these elements I will get
!         at least 0, at most nb and, if there is just part of a block
!         for me, I'll get rest-myrow*nb elements more
!         Thus the elements I get in addition to the (jglob-1)/(nb*nprow)*nb
!         each processor gets is
!         min(max(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
          irecvcnts(jprocrow)=irecvcnts(jprocrow)+&
             (jglob-1)/(nb*nprow)*nb + &
             min(max(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
       endif
    enddo  ! jglob=k,min(k+nb-1,n)
 enddo ! k=mycol*nb+1,n,npcol*nb
!  Compute length of new receivebuffer
!  and start addresses in the receivebuffer
 istartrecvbuf(0)=0
 iinsendbuf(0)=0
 do jprocrow=0,nprow-1
    if (jprocrow /= myrow) then ! I will get and send in this alltoall
       lenrecvbuf2=lenrecvbuf2+irecvcnts(jprocrow)
       if (jprocrow < nprow-1) then
          istartrecvbuf(jprocrow+1)=istartrecvbuf(jprocrow)+&
             irecvcnts(jprocrow)
          iinsendbuf(jprocrow+1)=iinsendbuf(jprocrow)+ &
             isendcnts(jprocrow)
       endif
    else  ! this is my row, i won't get anything but have to count
!              start in receivebuffer
       if (jprocrow < nprow-1) then
          istartrecvbuf(jprocrow+1)= istartrecvbuf(jprocrow)
          iinsendbuf(jprocrow+1)=iinsendbuf(jprocrow)
       endif
    endif
 enddo
 If (lenrecvbuf2 > lenrecvbuf) then
    Deallocate(recvbuffer) ! Old one no longer needed if it is too small
    Allocate(recvbuffer(lenrecvbuf2))
 Endif
#ifdef CPP_INVERSION
 Call MPI_ALLTOALLV(sendbuffer,isendcnts,iinsendbuf,&
       CPP_MPI_TYP_REAL,recvbuffer,irecvcnts,istartrecvbuf,&
       CPP_MPI_TYP_REAL,icommcol(mycol),ierr)
#else
 Call MPI_ALLTOALLV(sendbuffer,isendcnts,iinsendbuf,&
       CPP_MPI_TYP_COMPLEX,recvbuffer,irecvcnts,istartrecvbuf,&
       CPP_MPI_TYP_COMPLEX,icommcol(mycol),ierr)
#endif
 deallocate(sendbuffer)
 istartrecvbuf=istartrecvbuf+1
!
!  After second alltoall store elements received to Asca
!  Go again through colums I have to get
!
 do k=mycol*nb+1,n,npcol*nb ! k is the first element in a column block
    do jglob=k,min(k+nb-1,n) ! global column index
!      the processor originally having this column is 
!      processor modulo(jglob-1,np)
       If (modulo(jglob-1,np)/npcol /= myrow) then ! I received it in second alltoallv
!         modulo(jglob-1,np) is the number of the processor who had
!         this column in the beginning
!         This processor is in my processor row, if its number divided by
!         npcol is the same as my number divided by npcol, and this is
!         myrow 
!         In that case nothing has to be done, I already treated this
!         column before second alltoall
          jprocrow=modulo(jglob-1,np)/npcol
!         The number of elements I have to receive for this column is:
!         The number of complete "sweeps" of nb*nprow elements of this
!         column are (jglob-1)/(nb*nprow), that is, each processor in my
!         processor column gets at least (jglob-1)/(nb*nprow)*nb elements
!         of column jglob.
!         The rest is jglob-(jglob-1)/(nb*nprow)*nb*nprow elements.
!         Out of these elements I will get
!         at least 0, at most nb and, if there is just part of a block
!         for me, I'll get rest-myrow*nb elements more
!         Thus the elements I get in addition to the (jglob-1)/(nb*nprow)*nb
!         each processor gets is
!         min(max(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
          jloc=jglob-(jglob-1)/(nb*npcol)*nb*(npcol-1)-nb*mycol
          Do i=1,(jglob-1)/(nb*nprow)*nb + &
             min(max(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
             Asca(i,jloc)=recvbuffer(istartrecvbuf(jprocrow)+i-1)
          Enddo
          istartrecvbuf(jprocrow) = istartrecvbuf(jprocrow)+&
             (jglob-1)/(nb*nprow)*nb + &
             min(max(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
       endif ! End if I got it in second alltoallv
    enddo  ! jglob=k,min(k+nb-1,n)
 enddo ! k=mycol*nb+1,n,npcol*nb
!
! Free the row and column communicators
!
 do k=0,nprow-1
    if (icommrow(k) /= MPI_COMM_NULL) then
       call MPI_COMM_FREE(icommrow(k),ierr)
    endif
    call MPI_GROUP_FREE(irowgroup(k),ierr)
 enddo
 do k=0,npcol-1
    if (icommcol(k) /= MPI_COMM_NULL) then
       call MPI_COMM_FREE(icommcol(k),ierr)
    endif
    call MPI_GROUP_FREE(icolgroup(k),ierr)
 enddo
!
!  Deallocate the arrays to store the indices
!
 deallocate(irecvcnts) 
 deallocate(istartrecvbuf) 
 deallocate(iinsendbuf) 
 deallocate(indx) 
 deallocate(isendcnts)    
 deallocate(iwork)
 deallocate(irowgroup)
 deallocate(icommrow)
 deallocate(icolgroup)
 deallocate(icommcol)
return
end subroutine subredist1
