MODULE m_fermie
  USE m_juDFT
  !-----------------------------------------------------------------------
  !     determines the fermi energy by
  !            gaussian-integration method                          c.l.fu
  !            triangular method (or tetrahedrons)
  !            or fermi-function                                    p.kurz
  !----------------------------------------------------------------------
CONTAINS
  SUBROUTINE fermie(eig_id, mpi,kpts,obsolete,&
       input, noco,e_min,jij,cell,results)

    !---------------------------------------------------f--------------------
    !
    !     a fist (T=0) approximation to the fermi-energy is determined
    !     by:
    !           zelec = sum { spindg * we }
    !                       e=<ef
    !
    !     TREE STRUCTURE: fermie----sort
    !                           ----fergwt
    !                           ----fertri---triang
    !                                     |--dosint
    !                                     |--dosef -- trisrt
    !                                     +--doswt
    !                           ----ferhis---ef_newton
    !
    !-----------------------------------------------------------------------


    USE m_eig66_io, ONLY : read_eig
#if defined(CPP_MPI)&&defined(CPP_NEVER)
    USE m_mpi_col_eigJ
#endif
    USE m_sort
    USE m_fertri
    USE m_ferhis
    USE m_fergwt
    USE m_types
    USE m_xmlOutput
    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_mpi),INTENT(IN)   :: mpi
    TYPE(t_obsolete),INTENT(IN)   :: obsolete
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_jij),INTENT(IN)   :: jij
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_kpts),INTENT(IN)   :: kpts
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id
    REAL,INTENT(IN)      :: e_min
    !     ..
    !     .. Array Arguments ..
    !REAL,    INTENT (OUT):: w(:,:,:) !(dimension%neigd,kpts%nkpt,dimension%jspd)
    !     ..
    !     .. Local Scalars ..
    REAL del  ,spindg,ssc ,ws,zc,tkb_1,weight
    INTEGER i,idummy,j,jsp,k,l,n,nbands,nstef,nv,nmat,nspins
    INTEGER n_help
    !     ..
    !     .. Local Arrays ..
    !
    INTEGER, ALLOCATABLE :: idxeig(:),idxjsp(:),idxkpt(:),INDEX(:)
    REAL,    ALLOCATABLE :: e(:),eig(:,:,:),we(:)
    INTEGER ne(kpts%nkpt,SIZE(results%w_iks,3))
    CHARACTER(LEN=20)    :: attributes(5)

    !--- J constants
    !--- J constants

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER, PARAMETER :: comm = MPI_COMM_SELF
    INTEGER*4 :: nv_mpi(2),idum1d(0),idum2d(0,0)
#endif

    !     ..
    !     ..
    !***********************************************************************
    !                  ABBREVIATIONS
    !
    !     eig        : array of eigenvalues 
    !     wtkpt      : list of the weights of each k-point (from inp-file)
    !     e          : linear list of the eigenvalues 
    !     we         : list of weights of the eigenvalues in e
    !     zelec      : number of electrons 
    !     spindg     : spindegeneracy (2 in nonmagnetic calculations)
    !     seigv      : weighted sum of the occupied valence eigenvalues
    !     seigsc     : weighted sum of the semi-core eigenvalues
    !     seigscv    : sum of seigv and seigsc
    !     ts         : entropy contribution to the free energy
    !
    !***********************************************************************
    !     .. Data statements ..
    DATA del/1.0e-6/
    !     ..
    n=SIZE(results%w_iks) !size of list of all eigenvalues
    ALLOCATE (idxeig(n),idxjsp(n),idxkpt(n),INDEX(n),e(n),we(n) )
    ALLOCATE (eig(SIZE(results%w_iks,1),SIZE(results%w_iks,2),SIZE(results%w_iks,3)))

    ! initiliaze e
    e = 0
    !
    IF (jij%l_J) THEN
#if defined(CPP_MPI)&&defined(CPP_NEVER)
       CALL mpi_col_eigJ(mpi%mpi_comm,mpi%irank,mpi%isize,kpts%nkpt,SIZE(results%w_iks,1),kpts%nkpt,&
            &                       jij%nkpt_l,jij%eig_l,&
            &                    kpts%bk,kpts%wtkpt,ne(1,1),eig)
       IF (mpi%irank.NE.0) THEN
          DEALLOCATE( idxeig,idxjsp,idxkpt,index,e,eig,we )
          RETURN
       ENDIF
#endif
    ENDIF


    IF ( mpi%irank == 0 ) WRITE (6,FMT=8000)
8000 FORMAT (/,/,1x,'fermi energy and band-weighting factors:')
    !
    !---> READ IN EIGENVALUES
    !
    spindg = 2.0/REAL(input%jspins)
    n = 0
    results%seigsc = 0.0
    ssc = 0.0
    n_help = 0
    !
    !---> pk non-collinear
    IF (noco%l_noco) THEN
       nspins = 1
    ELSE
       nspins = input%jspins
    ENDIF
    !---> pk non-collinear
    !
    IF (mpi%irank == 0) CALL openXMLElementNoAttributes('eigenvalues')
    DO jsp = 1,nspins
       DO  k = 1,kpts%nkpt
          CALL read_eig(eig_id,k,jsp,neig=ne(k,jsp),eig=eig(:,k,jsp))
          IF ( mpi%irank == 0 ) THEN
             WRITE (6,'(a2,3f10.5,f12.6)') 'at',kpts%bk(:,k),kpts%wtkpt(k)
             WRITE (6,'(i5,a14)') ne(k,jsp),' eigenvalues :' 
             WRITE (6,'(8f12.6)') (eig(i,k,jsp),i=1,ne(k,jsp))
             attributes = ''
             WRITE(attributes(1),'(i0)') jsp
             WRITE(attributes(2),'(i0)') k
             WRITE(attributes(3),'(f15.8)') kpts%bk(1,k)
             WRITE(attributes(4),'(f15.8)') kpts%bk(2,k)
             WRITE(attributes(5),'(f15.8)') kpts%bk(3,k)
             CALL writeXMLElementPoly('eigenvaluesAt',(/'spin','ikpt','k_x ','k_y ','k_z '/),attributes,eig(1:ne(k,jsp),k,jsp))
          END IF
          nv= -1
          !
          !--->          STORE EIGENVALUES AND WEIGHTS IN A LINEAR LIST. AND MEMORIZE 
          !--->          CONECTION TO THE ORIGINAL ARRAYS
          !
          DO  j = 1,ne(k,jsp)
             e(n+j) = eig(j,k,jsp)
             we(n+j) = kpts%wtkpt(k)
             idxeig(n+j) = j+n_help
             idxkpt(n+j) = k
             idxjsp(n+j) = jsp
 	       END DO
          !--->          COUNT THE NUMBER OF EIGENVALUES
          n = n + ne(k,jsp)
       END DO
    END DO
    IF (mpi%irank == 0) CALL closeXMLElement('eigenvalues')

    CALL sort(n,e,index)

    !     Check if no deep eigenvalue is found
    IF (e_min-MINVAL(e(1:n))>1.0) THEN
       WRITE(6,*) 'WARNING: Too low eigenvalue detected:'
       WRITE(6,*) 'min E=', MINVAL(e(1:n)),' min(enpara)=',&
            &             e_min
       CALL juDFT_warn("Too low eigenvalue detected",calledby="fermi" &
            &     ,hint ="If the lowest eigenvalue is more than 1Htr below "//&
            &     "the lowest energy parameter, you probably have picked up"//&
            &     " a ghoststate")
    END IF
    !
    !---> DETERMINE EF BY SUMMING WEIGHTS
    !
    weight = input%zelec/spindg
    results%seigv = 0.0e0
    ws = 0.0e0
    l = 0
    DO WHILE ((ws+del).LT.weight)
       l = l + 1
       IF (l.GT.n) THEN
          IF ( mpi%irank == 0 ) THEN
             WRITE (16,FMT=8010) n,ws,weight
             WRITE (6,FMT=8010) n,ws,weight
          END IF
          CALL juDFT_error("fermi",calledby="fermie")
8010      FORMAT (/,10x,'error: not enough wavefunctions.',i10,&
               &             2d20.10)
       END IF
       ws = ws + we(INDEX(l))
       results%seigv = results%seigv + e(INDEX(l))*we(INDEX(l))*spindg
       !         WRITE (6,FMT='(2f10.7)') e(index(l)),we(index(l))
    END DO
    results%ef = e(INDEX(l))
    nstef = l
    zc = input%zelec
    IF ( mpi%irank == 0 ) WRITE (6,FMT=8020) results%ef,nstef,results%seigv,ws,results%seigsc,ssc

    !+po
    results%ts = 0.0
    !-po
    results%w_iks = 0.0
    IF (input%gauss) THEN
       CALL fergwt(kpts,input,mpi,ne, eig,results)
    ELSE IF (input%tria) THEN
       CALL fertri(input,kpts,mpi%irank, ne,kpts%nkpt,nspins,zc,eig,kpts%bk,spindg,&
            results%ef,results%seigv,results%w_iks)
    ELSE
       nspins = input%jspins
       IF (noco%l_noco) nspins = 1
       tkb_1 = input%tkb
       CALL ferhis(input,kpts,mpi,results,index,idxeig,idxkpt,idxjsp, n,&
            nstef,ws,spindg,weight,e,ne,we, noco,jij,cell)
    END IF
    !     7.12.95 r.pentcheva seigscv must be calculated outside if (gauss)
    results%seigscv = results%seigsc + results%seigv
    !
    DEALLOCATE ( idxeig,idxjsp,idxkpt,index,e,eig,we )

    attributes = ''
    WRITE(attributes(1),'(f20.10)') results%ef
    WRITE(attributes(2),'(a)') 'Htr'
    IF (mpi%irank.EQ.0) CALL writeXMLElement('FermiEnergy',(/'value','units'/),attributes(1:2))

    RETURN
8020 FORMAT (/,'FERMIE:',/,&
         &       10x,'first approx. to ef    (T=0)  :',f10.6,' htr',&
         &       '   (energy of the highest occ. eigenvalue)',/,&
         &       10x,'number of occ. states  (T=0)  :',i10,/,&
         &       10x,'first approx. to seigv (T=0)  :',f10.6,' htr',/,&
         &       10x,'sum of weights of occ. states :',f10.6,/,&
         &       10x,'sum of semicore eigenvalues   :',f10.6,' htr',/,&
         &       10x,'sum of semicore charge        :',f10.6,' e',/)
  END SUBROUTINE fermie
     END MODULE m_fermie
