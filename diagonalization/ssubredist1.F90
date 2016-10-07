!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!DEC$ FREEFORM
SUBROUTINE subredist1(n,lda,SUB_COMM,nprow,npcol,iam,ierr,nb,achi_r,asca_r,achi_c,asca_c)
  USE m_juDFT
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
  REAL,OPTIONAL     :: achi_r(:) ! Input, matrix in chani-distribution
  REAL,OPTIONAL     :: asca_r(lda,*) ! Output, matrix in ScaLAPACK distribution
  COMPLEX,OPTIONAL  :: achi_c(*) ! Input, matrix in chani-distribution
  COMPLEX,OPTIONAL  :: asca_c(lda,*) ! Output, matrix in ScaLAPACK distribution
  ! Matrix might be real or complex 

  INTEGER  :: n,lda   ! Global matrix size, local leading dimension of asca
  INTEGER  :: SUB_COMM,nprow,npcol,iam ! Communicator, 
  !                number of processor rows and columns in SUB_COMM, my rank
  INTEGER  :: nb   ! blocking size, will be computed later
  INTEGER  :: ierr  ! Error parameter to report problems
  !
  ! Local Variables
  INTEGER :: np       ! number of procs, nprow*npcol=np
  INTEGER :: me, myrow, mycol ! my processor coordinates
  INTEGER :: j,i,k,l, kk           ! counters
  INTEGER :: mynumcols             ! number of columns processor iam owns in the beginning
  INTEGER :: jproccol, jprocrow    ! processor column/row number counter variable
  INTEGER :: isndproc, him         ! processor in my column sending to me
  INTEGER :: hisnumcols            ! number of matrix columns proc isndproc owns
  INTEGER :: iworldgroup           ! mpi group index for group to SUB_COMM
  INTEGER :: myrowssca,mycolssca   ! number of rows/cols processor iam gets in Asca
  INTEGER, ALLOCATABLE :: irowgroup(:), icolgroup(:) ! Subgroups rows, columns
  INTEGER, ALLOCATABLE :: icommrow(:), icommcol(:) ! Communicators for Subgroups
  INTEGER, ALLOCATABLE :: iranks(:)  ! processor ranks to create subgroups
  INTEGER, ALLOCATABLE :: glcolnum(:,:)  ! Global column number of the columns I get
  INTEGER :: jstart, jend, npfac   ! see comments where used
  INTEGER, ALLOCATABLE :: iwork(:,:), indx(:), isendcnts(:)  ! see comments where used
  INTEGER, ALLOCATABLE :: irecvcnts(:)  ! see comments where used
  INTEGER, ALLOCATABLE :: irecvcol(:)  ! see comments where used
  !Integer, allocatable :: ihavecolnum(:)  ! Which columns do I have after first alltoall
  INTEGER, ALLOCATABLE :: istartrecvbuf(:)  ! Start of cols from proc i in recvbuf
  INTEGER, ALLOCATABLE :: iinsendbuf(:)  ! Now in sendbuf for proc i
  INTEGER :: nsq      ! Integer sqrt of np, to compute grid
  INTEGER :: lensendbuf, lenrecvbuf ! Length of complete send- and recvbuf, first send
  INTEGER :: lensendbuf2, lenrecvbuf2 ! Length of complete send- and recvbuf, second send
  INTEGER :: jglob, jloc    ! Global and local column index
  INTEGER :: lglob, lloc    ! Global and local row index
  !Integer :: icntcol    ! Number of columns I own after first alltoall
  INTEGER :: indxrecvbuf    ! Where I am in receivebuffer of first alltoall
  INTEGER :: ierr1    ! For all the MPI calls
  INTEGER :: subsize,mysubrank    ! Test of comm_create
  ! Arrays to redistribute
  !
  REAL,ALLOCATABLE :: sendbuffer_r(:)
  REAL,ALLOCATABLE :: recvbuffer_r(:)
  COMPLEX,ALLOCATABLE :: sendbuffer_c(:)
  COMPLEX,ALLOCATABLE :: recvbuffer_c(:)

  LOGICAL :: l_real
  l_real=PRESENT(achi_r)
  !
  ! Now do the computation
  !
  np=nprow*npcol   ! number of processors
  mynumcols=n/np   ! my number of columns in the chani-distribution
  IF (iam < n-(n/np)*np) mynumcols=mynumcols+1
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
  CALL MPI_COMM_GROUP(SUB_COMM,iworldgroup, Ierr)
  IF (ierr /= 0) THEN
     WRITE(6,*) 'MPI_COMM_GROUP failed, ierr=',ierr
     !    CALL CPP_flush(6)
     CALL juDFT_error("MPI_COMM_GROUP failed")
  ENDIF
  !
  ! Allocate arrays to build subgroups for rows and columns
  ! and create subgroups consisting of processor rows and columns
  !
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
     IF (ierr /= 0) THEN
        WRITE(6,*) 'MPI_GROUP_INCL(',k,') row failed, ierr=',ierr
        !       CALL CPP_flush(6)
        CALL juDFT_error("MPI_COMM_INCL failed")
     ENDIF
     CALL MPI_COMM_CREATE(SUB_COMM,irowgroup(k),&
          icommrow(k),ierr)
     IF (ierr /= 0) THEN
        WRITE(6,*) 'MPI_COMM_CREATE(',k,') row failed, ierr=',ierr
        !       CALL CPP_flush(6)
        CALL juDFT_error("MPI_COMM_CREATE failed")
     ENDIF
  ENDDO
  DO k=0,npcol-1  ! Create groups of one processor column each
     DO i=1,nprow
        iranks(i)=k+(i-1)*npcol
     ENDDO
     CALL MPI_GROUP_INCL(iworldgroup,nprow,iranks,icolgroup(k),&
          ierr)
     IF (ierr /= 0) THEN
        WRITE(6,*) 'MPI_GROUP_INCL(',k,') col failed, ierr=',ierr
        !CALL CPP_flush(6)
        CALL juDFT_error("MPI_COMM_INCL failed")
     ENDIF
     CALL MPI_COMM_CREATE(SUB_COMM,icolgroup(k),&
          icommcol(k),ierr)
     IF (ierr /= 0) THEN
        WRITE(6,*) 'MPI_COMM_CREATE(',k,') col failed, ierr=',ierr
        !CALL CPP_flush(6)
        CALL juDFT_error("MPI_COMM_CREATE failed")
     ENDIF
  ENDDO
  DEALLOCATE(iranks)
  !
  !
  ! Allocate the arrays to store the indices
  !
  ALLOCATE(indx(0:npcol-1))
  ALLOCATE(isendcnts(0:MAX(npcol-1,nprow-1)))
  ALLOCATE(iinsendbuf(0:MAX(npcol-1,nprow-1)))
  ! indx(jproccol) will indicate how many columns of processor column
  ! jproccol are already found on this processor
  indx=0  ! I have not found any column of any other processor yet
  isendcnts=0  ! The length is 0 thus
  ALLOCATE(iwork(1:2*mynumcols,0:npcol-1))
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
  DO j=1,mynumcols  ! Find out who my columns belong to
     jstart=jend+1  ! Start point in achi 
     jend=jend+me+npfac*np ! End point in achi
     jglob=me+npfac*np  ! global column number
     npfac=npfac+1  ! next column is np elements longer
     jproccol=MODULO(jglob-1,npcol*nb)/nb ! Belongs to proc column
     iwork(2*indx(jproccol)+1,jproccol)=jstart
     iwork(2*indx(jproccol)+2,jproccol)=jend
     isendcnts(jproccol)=isendcnts(jproccol)+jend-jstart+1
     !   Send isendcnts(jproccol) elements to (myrow,jproccol)
     indx(jproccol)=indx(jproccol)+1 ! Send indx(jproccol) columns to him
  ENDDO  ! End j-loop over my columns
  ! Count the number of columns and elements I have to send to which processor column
  ! Here I only want to compute lensendbuf, thus I don't need to store
  ! anything to sendbuffer because sendbuffer is not yet allocated!!!
  lensendbuf=0 ! Length of sendbuffer
  DO j=0,npcol-1  
     IF (j /= mycol) lensendbuf=lensendbuf+isendcnts(j)
  ENDDO
  IF (l_real) THEN
     ALLOCATE(sendbuffer_r(lensendbuf))
  ELSE
     ALLOCATE(sendbuffer_c(lensendbuf))
  ENDIF

  !
  !  Now I have to put the elements to be sent to sendbuffer !!
  !
  lensendbuf=0
  iinsendbuf(0)=0 ! Starting addresses in sendbuffer (See MPI_ALLTOALLV)
  DO j=0,npcol-1  
     IF (indx(j)==0) THEN
        IF (isendcnts(j) /= 0) WRITE(6,*) iam, &
             ' Fehler, isendcnts nicht 0',isendcnts(j)
        IF (j < npcol-1) iinsendbuf(j+1)=iinsendbuf(j)
     ELSE
        IF (j /= mycol) THEN
           IF (l_real) THEN
              DO i=1,indx(j)
                 DO kk=iwork(2*i-1,j),iwork(2*i,j)
                    sendbuffer_r(lensendbuf+1+kk-iwork(2*i-1,j))=&
                         Achi_r(kk)
                 ENDDO
                 lensendbuf=lensendbuf+iwork(2*i,j)-iwork(2*i-1,j)+1
              ENDDO
           ELSE
              DO i=1,indx(j)
                 DO kk=iwork(2*i-1,j),iwork(2*i,j)
                    sendbuffer_c(lensendbuf+1+kk-iwork(2*i-1,j))=&
                         Achi_c(kk)
                 ENDDO
                 lensendbuf=lensendbuf+iwork(2*i,j)-iwork(2*i-1,j)+1
              ENDDO
           ENDIF
        ENDIF
        IF (j == mycol) isendcnts(j)=0 ! must be 0 for mpi_alltoallv
        IF (j < npcol-1)  iinsendbuf(j+1)=iinsendbuf(j)+isendcnts(j)
     ENDIF
  ENDDO
  ALLOCATE(irecvcnts(0:MAX(npcol-1,nprow-1))) ! Number of elements I will get
  ALLOCATE(irecvcol(0:npcol-1)) ! Number of columns I will get
  ALLOCATE(glcolnum(n/np+1,0:npcol-1)) ! Global column numbers
  ALLOCATE(istartrecvbuf(0:MAX(npcol-1,nprow-1))) ! Start in receivebuffer
  irecvcnts=0
  irecvcol=0
  glcolnum=0
  lenrecvbuf=0
  istartrecvbuf(0)=0
  DO k=0,npcol-1 ! Now look what processors in my row send to me
     isndproc=myrow*npcol+k ! Number of processor (myrow,k)
     IF (k /= mycol) THEN  ! All except me are interesting
        him=isndproc+1  ! make life easier
        hisnumcols=n/np   ! his number of columns in the chani-distribution
        IF (isndproc < n-(n/np)*np) hisnumcols=hisnumcols+1
        jend=0  ! jstart will be jend+1
        npfac=0 ! his first elements belong to a column shorter or equal np
        DO j=1,hisnumcols  ! Look through his columns
           jstart=jend+1 ! Starting in his part of Achi
           jend=jend+him+npfac*np ! Ending in his part of Achi
           jglob=him+npfac*np ! Global column number
           npfac=npfac+1  ! next column is np elements longer
           jproccol=MODULO(jglob-1,npcol*nb)/nb ! Will belong to proc column
           IF (jproccol==mycol) THEN ! I will get this column
              irecvcol(k)=irecvcol(k)+1 ! I get one more column from (myrow,k)
              irecvcnts(k)=irecvcnts(k)+jend-jstart+1
              glcolnum(irecvcol(k),k)=jglob ! It is global column number
           ENDIF  ! I get this column
        ENDDO       ! j=1,hisnumcols
        lenrecvbuf=lenrecvbuf+irecvcnts(k)
        IF (k < npcol-1) istartrecvbuf(k+1) = &
             istartrecvbuf(k)+irecvcnts(k)  
        !         what I got from k starts at istartrecvbuf(k)+1 in the recvbuf
     ELSE ! it is me, start in recvbuf for the next is same as for this 
        IF (k < npcol-1) istartrecvbuf(k+1)= istartrecvbuf(k)
     ENDIF       ! If it is not me
  ENDDO       ! Loop over procs in my row
  IF (l_real) THEN
     ALLOCATE(recvbuffer_r(lenrecvbuf))
     CALL MPI_ALLTOALLV(sendbuffer_r,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_REAL,recvbuffer_r,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_REAL,icommrow(myrow),ierr)
  ELSE
     ALLOCATE(recvbuffer_c(lenrecvbuf))
     CALL MPI_ALLTOALLV(sendbuffer_c,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_COMPLEX,recvbuffer_c,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_COMPLEX,icommrow(myrow),ierr)
  ENDIF
  !
  !  Now look to Asca to put the elements of Achi or receivebuffer to
  !  See later, where the number of elements of column jglob this
  !  processor will get is computed. Jglob there plays the same role as
  !  n here.
  !
  myrowssca=(n-1)/(nb*nprow)*nb+ &
       MIN(MAX(n-(n-1)/(nb*nprow)*nb*nprow-nb*myrow,0),nb)
  !  My number of rows in ScaLAPACK distribution
  IF (myrowssca > lda) &
       WRITE(6,*) 'Redist1: Wrong dimension of Asca, lda=',lda,  &
       ' myrowssca=',myrowssca
  ! My number of columns in ScaLAPACK distribution
  mycolssca=(n-1)/(nb*npcol)*nb+ &
       MIN(MAX(n-(n-1)/(nb*npcol)*nb*npcol-nb*mycol,0),nb)
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
  DO k=1,indx(mycol)    ! I had indx(mycol) columns for my processor column
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
     DO j=1,jglob,nb
        !      Go through column jglob by row blocks
        !      The row block treated now belongs to processor row
        !      mod(j/nb,nprow)
        IF (MOD(j/nb,nprow)==myrow) THEN ! Belonging to my processor row 
           !                                         = I get it, store it to Asca
           IF (l_real) THEN
              DO lglob=j,MIN(jglob,j+nb-1) ! This row block, no more than jglob
                 Asca_r(lloc,jloc)=Achi_r(iwork(2*k-1,mycol)+lglob-1)
                 lloc=lloc+1
              ENDDO
           ELSE
              DO lglob=j,MIN(jglob,j+nb-1) ! This row block, no more than jglob
                 Asca_c(lloc,jloc)=Achi_c(iwork(2*k-1,mycol)+lglob-1)
                 lloc=lloc+1
              ENDDO
           END IF
        ELSE ! send to processor row mod(j/nb,nprow)
           isendcnts(MOD(j/nb,nprow))=isendcnts(MOD(j/nb,nprow))&
                + MIN(nb,jglob-j+1) ! Compute only number of elements
        ENDIF
     ENDDO
  ENDDO
  !
  !  Now processing the received columns, processor by processor
  !
  ! icntcol=indx(mycol)
  indxrecvbuf=0  ! Go through receivebuffer element by element
  DO l=0,npcol-1
     DO k=1,irecvcol(l) !  If l=mycol irecvcol=0, nothing will be done
        jglob=glcolnum(k,l) ! Received column is global column number jglob
        jloc=jglob-(jglob-1)/(nb*npcol)*nb*(npcol-1)-nb*mycol
        !      It will be local column number jloc of Asca
        !       icntcol = icntcol+1 ! Number of columns I have now
        !       ihavecolnum(icntcol)=jglob ! Global indices of columns I have
        lloc=1 ! Now go through column by row blocks
        DO j=1,jglob,nb
           !         Go through column jglob by row blocks
           !         The row block treated now belongs to processor row
           !         mod(j/nb,nprow)
           IF (MOD(j/nb,nprow)==myrow) THEN  ! my row block, store it to Asca
              IF (l_real) THEN
                 DO lglob=j,MIN(jglob,j+nb-1) ! This row block, no more than jglob
                    indxrecvbuf=indxrecvbuf+1
                    Asca_r(lloc,jloc)=recvbuffer_r(indxrecvbuf)
                    lloc=lloc+1
                 ENDDO
              ELSE
                 DO lglob=j,MIN(jglob,j+nb-1) ! This row block, no more than jglob
                    indxrecvbuf=indxrecvbuf+1
                    Asca_c(lloc,jloc)=recvbuffer_c(indxrecvbuf)
                    lloc=lloc+1
                 ENDDO
              ENDIF
           ELSE  ! not my row block send to processor row mod(j/nb,nprow)
              !                 min(nb,jglob-j+1) more elements
              isendcnts(MOD(j/nb,nprow))=   &
                   isendcnts(MOD(j/nb,nprow))+ MIN(nb,jglob-j+1)
              indxrecvbuf=indxrecvbuf+MIN(nb,jglob-j+1)
           ENDIF ! who's row block
        ENDDO ! End going through column jglob
     ENDDO ! End going through columns received from processor myrow,k
  ENDDO ! End going through all processors in my row
  ! Deallocate some workspace no longer needed
  DEALLOCATE(glcolnum)
  DEALLOCATE(irecvcol)
  !  deallocate(ihavecolnum)
  !
  iinsendbuf(0)=1
  DO k=0,nprow-1
     IF (k /= myrow) THEN ! I will send isendcnts(k )elements to processor
        !                          k,mycol
        lensendbuf2=lensendbuf2+isendcnts(k)
        IF (k < nprow-1) iinsendbuf(k+1)=iinsendbuf(k)+isendcnts(k)
     ELSE  ! It is my row, don't have to send
        IF (k < nprow-1) iinsendbuf(k+1)=iinsendbuf(k)
     ENDIF
  ENDDO
  IF (lensendbuf2 > lensendbuf) THEN
     IF (l_real) THEN
        DEALLOCATE(sendbuffer_r)
        ALLOCATE(sendbuffer_r(lensendbuf2))
     ELSE
        DEALLOCATE(sendbuffer_c)
        ALLOCATE(sendbuffer_c(lensendbuf2))
     ENDIF
  ENDIF
  !
  !  Now look through all columns I need and see whether I have them or get them
  !
  irecvcnts=0 ! Count for next receive
  lenrecvbuf2=0  ! New receivebuffer 
  jstart=0    ! No column from Achi stored to sendbuf
  DO k=mycol*nb+1,n,npcol*nb ! k is the first element in a column block
     DO jglob=k,MIN(k+nb-1,n) ! global column index
        !      the processor originally having this column is 
        !      processor modulo(jglob-1,np)
        IF (MODULO(jglob-1,np)==iam) THEN ! It was my column in the beginning
           jstart=jstart+1
           DO j=1,jglob,nb ! Row block starting with j
              !            Go through column jglob by row blocks
              !            The row block treated now belongs to processor row
              !            mod(j/nb,nprow)
              IF (MOD(j/nb,nprow) /= myrow) THEN  ! not my row block
                 !               Now find out where it is in Achi
                 IF (l_real) THEN
                    DO l=j,MIN(j+nb-1,jglob)
                       sendbuffer_r(iinsendbuf(MOD(j/nb,nprow)))=&
                            Achi_r(iwork(2*jstart-1,mycol)+l-1)
                       iinsendbuf(MOD(j/nb,nprow))= &
                            iinsendbuf(MOD(j/nb,nprow))+1
                    ENDDO
                 ELSE
                    DO l=j,MIN(j+nb-1,jglob)
                       sendbuffer_c(iinsendbuf(MOD(j/nb,nprow)))=&
                            Achi_c(iwork(2*jstart-1,mycol)+l-1)
                       iinsendbuf(MOD(j/nb,nprow))= &
                            iinsendbuf(MOD(j/nb,nprow))+1
                    ENDDO
                 ENDIF
              ENDIF ! who's row block
           ENDDO ! End going through column jglob
        ELSE IF (MODULO(jglob-1,np)/npcol==myrow) THEN
           jproccol=MODULO(MODULO(jglob-1,np),npcol)
           !         modulo(jglob-1,np) is the number of the processor who had
           !         this column in the beginning
           !         This processor is in my processor row, if its number divided by
           !         npcol is the same as my number divided by npcol, and this is
           !         myrow
           !         processor number divided by npcol = processor row number
           !         processor number modulo npcol = processor column number
           !         It is in my receive buffer
           DO j=1,jglob,nb
              !            Go through column jglob by row blocks
              !            The row block treated now belongs to processor row
              !            mod(j/nb,nprow)
              IF (MOD(j/nb,nprow)==myrow) THEN  ! my row block
                 !               stored the elements from recvbuf to Asca
                 istartrecvbuf(jproccol)=istartrecvbuf(jproccol)+&
                      MIN(nb,jglob-j+1)
              ELSE  ! not my row block, store to sendbuffer for row mod(j/nb,nprow)
                 !               Now find out where it is in recvbuffer
                 DO l=j,MIN(j+nb-1,jglob)
                    IF (l_real) THEN
                       sendbuffer_r(iinsendbuf(MOD(j/nb,nprow)))=recvbuffer_r(istartrecvbuf(jproccol))
                    ELSE
                       sendbuffer_c(iinsendbuf(MOD(j/nb,nprow)))=recvbuffer_c(istartrecvbuf(jproccol))
                    END IF

                    iinsendbuf(MOD(j/nb,nprow))= &
                         iinsendbuf(MOD(j/nb,nprow))+1
                    istartrecvbuf(jproccol)= &
                         istartrecvbuf(jproccol)+1
                 ENDDO
              ENDIF ! who's row block
           ENDDO ! End going through column jglob
        ELSE  ! I get it in this alltoall
           jprocrow=MODULO(jglob-1,np)/npcol ! processor row where it is
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
                MIN(MAX(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
        ENDIF
     ENDDO  ! jglob=k,min(k+nb-1,n)
  ENDDO ! k=mycol*nb+1,n,npcol*nb
  !  Compute length of new receivebuffer
  !  and start addresses in the receivebuffer
  istartrecvbuf(0)=0
  iinsendbuf(0)=0
  DO jprocrow=0,nprow-1
     IF (jprocrow /= myrow) THEN ! I will get and send in this alltoall
        lenrecvbuf2=lenrecvbuf2+irecvcnts(jprocrow)
        IF (jprocrow < nprow-1) THEN
           istartrecvbuf(jprocrow+1)=istartrecvbuf(jprocrow)+&
                irecvcnts(jprocrow)
           iinsendbuf(jprocrow+1)=iinsendbuf(jprocrow)+ &
                isendcnts(jprocrow)
        ENDIF
     ELSE  ! this is my row, i won't get anything but have to count
        !              start in receivebuffer
        IF (jprocrow < nprow-1) THEN
           istartrecvbuf(jprocrow+1)= istartrecvbuf(jprocrow)
           iinsendbuf(jprocrow+1)=iinsendbuf(jprocrow)
        ENDIF
     ENDIF
  ENDDO
  IF (lenrecvbuf2 > lenrecvbuf) THEN
     IF (l_real) THEN
        DEALLOCATE(recvbuffer_r) ! Old one no longer needed if it is too small
        ALLOCATE(recvbuffer_r(lenrecvbuf2))
     ELSE
        DEALLOCATE(recvbuffer_c) ! Old one no longer needed if it is too small
        ALLOCATE(recvbuffer_c(lenrecvbuf2))
     ENDIF
  ENDIF
  IF (l_real) THEN
     CALL MPI_ALLTOALLV(sendbuffer_r,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_REAL,recvbuffer_r,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_REAL,icommcol(mycol),ierr)
  ELSE
     CALL MPI_ALLTOALLV(sendbuffer_c,isendcnts,iinsendbuf,&
          CPP_MPI_TYP_COMPLEX,recvbuffer_c,irecvcnts,istartrecvbuf,&
          CPP_MPI_TYP_COMPLEX,icommcol(mycol),ierr)
  ENDIF
  IF (l_real) THEN
     DEALLOCATE(sendbuffer_r)
  ELSE
     DEALLOCATE(sendbuffer_c)
  END IF

  istartrecvbuf=istartrecvbuf+1
  !
  !  After second alltoall store elements received to Asca
  !  Go again through colums I have to get
  !
  DO k=mycol*nb+1,n,npcol*nb ! k is the first element in a column block
     DO jglob=k,MIN(k+nb-1,n) ! global column index
        !      the processor originally having this column is 
        !      processor modulo(jglob-1,np)
        IF (MODULO(jglob-1,np)/npcol /= myrow) THEN ! I received it in second alltoallv
           !         modulo(jglob-1,np) is the number of the processor who had
           !         this column in the beginning
           !         This processor is in my processor row, if its number divided by
           !         npcol is the same as my number divided by npcol, and this is
           !         myrow 
           !         In that case nothing has to be done, I already treated this
           !         column before second alltoall
           jprocrow=MODULO(jglob-1,np)/npcol
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
           IF (l_real) THEN
              DO i=1,(jglob-1)/(nb*nprow)*nb + &
                   MIN(MAX(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
                 Asca_r(i,jloc)=recvbuffer_r(istartrecvbuf(jprocrow)+i-1)
              ENDDO
           ELSE
              DO i=1,(jglob-1)/(nb*nprow)*nb + &
                   MIN(MAX(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
                 Asca_c(i,jloc)=recvbuffer_c(istartrecvbuf(jprocrow)+i-1)
              ENDDO
           END IF
           istartrecvbuf(jprocrow) = istartrecvbuf(jprocrow)+&
                (jglob-1)/(nb*nprow)*nb + &
                MIN(MAX(0,jglob-(jglob-1)/(nb*nprow)*nb*nprow-myrow*nb),nb)
        ENDIF ! End if I got it in second alltoallv
     ENDDO  ! jglob=k,min(k+nb-1,n)
  ENDDO ! k=mycol*nb+1,n,npcol*nb
  !
  ! Free the row and column communicators
  !
  DO k=0,nprow-1
     IF (icommrow(k) /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(icommrow(k),ierr)
     ENDIF
     CALL MPI_GROUP_FREE(irowgroup(k),ierr)
  ENDDO
  DO k=0,npcol-1
     IF (icommcol(k) /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(icommcol(k),ierr)
     ENDIF
     CALL MPI_GROUP_FREE(icolgroup(k),ierr)
  ENDDO
  !
  !  Deallocate the arrays to store the indices
  !
  DEALLOCATE(irecvcnts) 
  DEALLOCATE(istartrecvbuf) 
  DEALLOCATE(iinsendbuf) 
  DEALLOCATE(indx) 
  DEALLOCATE(isendcnts)    
  DEALLOCATE(iwork)
  DEALLOCATE(irowgroup)
  DEALLOCATE(icommrow)
  DEALLOCATE(icolgroup)
  DEALLOCATE(icommcol)
  RETURN
END SUBROUTINE subredist1
