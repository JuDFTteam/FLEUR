MODULE m_fermie
  USE m_juDFT
#ifdef CPP_MPI 
   use mpi 
#endif 
  !-----------------------------------------------------------------------
  !     determines the fermi energy by
  !            gaussian-integration method                          c.l.fu
  !            triangular method (or tetrahedrons)
  !            or fermi-function                                    p.kurz
  !----------------------------------------------------------------------
CONTAINS
  SUBROUTINE fermie(eig_id,fmpi,kpts,input,noco,e_min,cell,results)

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

    USE m_types
    USE m_constants
    USE m_eig66_io, ONLY : read_eig,write_eig
    USE m_sort
    USE m_fertri
    USE m_ferhis
    USE m_fergwt
    USE m_fertetra
    USE m_xmlOutput

    IMPLICIT NONE

    TYPE(t_results), INTENT(INOUT) :: results
    TYPE(t_mpi),     INTENT(IN)    :: fmpi
    TYPE(t_input),   INTENT(IN)    :: input
    TYPE(t_noco),    INTENT(IN)    :: noco
    TYPE(t_cell),    INTENT(IN)    :: cell
    TYPE(t_kpts),    INTENT(IN)    :: kpts
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id
    REAL,INTENT(IN)      :: e_min
    !     ..
    !     .. Array Arguments ..
    !REAL,    INTENT (OUT):: w(:,:,:) !(input%neig,kpts%nkpt,dimension%jspd)
    !     ..
    !     .. Local Scalars ..
    REAL del  ,spindg,ssc ,ws,zc,weight,efermi,seigv
    INTEGER i,idummy,j,jsp,k,l,n,nbands,nstef,nv,nmat,nspins
    INTEGER n_help,m_spins,mspin,sslice(2)
    !     ..
    !     .. Local Arrays ..
    !
    INTEGER :: idxeig(SIZE(results%w_iks)),idxjsp(SIZE(results%w_iks)),idxkpt(SIZE(results%w_iks)),INDEX(SIZE(results%w_iks))
    REAL    :: e(SIZE(results%w_iks)),we(SIZE(results%w_iks))
    CHARACTER(LEN=20)    :: attributes(5)

    !--- J constants
    !--- J constants

#ifdef CPP_MPI
    INTEGER, PARAMETER :: comm = MPI_COMM_SELF
    INTEGER*4 :: nv_mpi(2)
    INTEGER ierr
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
    !     ts         : entropy contribution to the free energy
    !
    !***********************************************************************
    !     .. Data statements ..
    DATA del/1.0e-6/

    ! initiliaze e
    e = 0

    IF ( fmpi%irank == 0 ) WRITE (oUnit,FMT=8000)
8000 FORMAT (/,/,1x,'fermi energy and band-weighting factors:')
    !
    !---> READ IN EIGENVALUES
    !
    spindg = 2.0/REAL(input%jspins)
    n = 0
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
    IF (fmpi%irank == 0. .and. .NOT.judft_was_argument("-minimalOutput")) CALL openXMLElementNoAttributes('eigenvalues')
    DO jsp = 1,nspins
       DO  k = 1,kpts%nkpt
          IF (fmpi%irank == 0) THEN
             if(input%eig66(1))CALL read_eig(eig_id,k,jsp,neig=results%neig(k,jsp),eig=results%eig(:,k,jsp))
             WRITE (oUnit,'(a2,3f10.5,a12,f12.6)') 'at',kpts%bk(:,k), " k-weight:", kpts%wtkpt(k)
             WRITE (oUnit,'(i5,a14)') results%neig(k,jsp),' eigenvalues :'
             WRITE (oUnit,'(8f12.6)') (results%eig(i,k,jsp),i=1,results%neig(k,jsp))
             IF(.NOT.judft_was_argument("-minimalOutput")) THEN
                attributes = ''
                WRITE(attributes(1),'(i0)') jsp
                WRITE(attributes(2),'(i0)') k
                WRITE(attributes(3),'(f15.8)') kpts%bk(1,k)
                WRITE(attributes(4),'(f15.8)') kpts%bk(2,k)
                WRITE(attributes(5),'(f15.8)') kpts%bk(3,k)
                CALL writeXMLElementPoly('eigenvaluesAt',(/'spin','ikpt','k_x ','k_y ','k_z '/),attributes,results%eig(1:results%neig(k,jsp),k,jsp))
             END IF
          END IF
#ifdef CPP_MPI
          CALL MPI_BARRIER(fmpi%mpi_comm,ierr)
#endif
       END DO
    ENDDO
    !finished reading of eigenvalues
    IF (fmpi%irank == 0 .and. .NOT.judft_was_argument("-minimalOutput")) CALL closeXMLElement('eigenvalues')
  IF (fmpi%irank == 0) THEN

    IF (ABS(input%fixed_moment)<1E-6) THEN
       !this is a standard calculation
       m_spins=1
    else
       !total moment is fixed
       m_spins=2
    END IF

    results%seigv = 0.0e0
    do mspin=1,m_spins
       IF (m_spins    == 1) THEN
          sslice = (/1,nspins/)
       ELSE
          sslice = (/mspin,mspin/)
          nspins = 1
       ENDIF
       n = 0
       DO jsp = sslice(1),sslice(2)
          !Generate a list of energies
          DO  k = 1,kpts%nkpt
             !
             !--->          STORE EIGENVALUES AND WEIGHTS IN A LINEAR LIST. AND MEMORIZE
             !--->          CONECTION TO THE ORIGINAL ARRAYS
             !
             DO  j = 1,results%neig(k,jsp)
                e(n+j) = results%eig(j,k,jsp)
                we(n+j) = kpts%wtkpt(k)
                idxeig(n+j) = j+n_help
                idxkpt(n+j) = k
                idxjsp(n+j) = jsp
             END DO
             !--->          COUNT THE NUMBER OF EIGENVALUES
             n = n + results%neig(k,jsp)
          END DO
       END DO

       CALL sort(index(:n),e)

       !     Check if no deep eigenvalue is found
       IF (e_min-MINVAL(e(1:n))>1.0) THEN
          WRITE(oUnit,*) 'WARNING: Too low eigenvalue detected:'
          WRITE(oUnit,*) 'min E=', MINVAL(e(1:n)),' min(enpara)=',e_min
          CALL juDFT_warn("Too low eigenvalue detected",calledby="fermi", &
                          hint ="If the lowest eigenvalue is more than 1Htr below "//&
                                "the lowest energy parameter, you probably have picked up"//&
                                " a ghoststate")
       END IF
       !
       !---> DETERMINE EF BY SUMMING WEIGHTS
       !
       weight = input%zelec/spindg
       seigv=0.0
       IF(m_spins /= 1) weight = weight/2.0  -(mspin-1.5)*input%fixed_moment
       ws = 0.0e0
       l = 0
       DO WHILE ((ws+del).LT.weight)
          l = l + 1
          IF (l.GT.n) THEN
             IF ( fmpi%irank == 0 ) THEN
                WRITE (oUnit,FMT=8010) n,ws,weight
             END IF
             CALL juDFT_error("Not enough wavefunctions",calledby="fermie")
8010         FORMAT (/,10x,'error: not enough wavefunctions.',i10,2d20.10)
          END IF
          ws = ws + we(INDEX(l))
          seigv =seigv + e(INDEX(l))*we(INDEX(l))*spindg
          !         WRITE (oUnit,FMT='(2f10.7)') e(index(l)),we(index(l))
       END DO
       results%ef = -100000.0
       IF(l.GT.0) THEN
          results%ef = e(INDEX(l))
       END IF
       nstef = l
       zc = input%zelec
       IF(m_spins /= 1) THEN
          zc = zc/2.0-(mspin-1.5)*input%fixed_moment
          idxjsp = 1 !assume single spin in following calculations
          IF (mspin == 1) THEN
             WRITE(oUnit,*) "Fixed total moment calculation"
             WRITE(oUnit,*) "Moment:",input%fixed_moment
             write(oUnit,*) "First Spin:"
          ELSE
             WRITE(oUnit,*) "Second Spin:"
          ENDIF
       ENDIF

       IF ( fmpi%irank == 0 ) WRITE (oUnit,FMT=8020) results%ef,nstef,seigv,ws,ssc

       !+po
       results%ts = 0.0
       !-po
       results%w_iks(:,:,sslice(1):sslice(2)) = 0.0
       results%bandgap = 0.0
       IF(input%bz_integration==BZINT_METHOD_HIST) THEN
          CALL ferhis(input,kpts,fmpi,index,idxeig,idxkpt,idxjsp,nspins, n,&
               nstef,ws,spindg,weight,e,results%neig(:,sslice(1):sslice(2)),we, noco,cell,results%ef,results%seigv,results%w_iks(:,:,sslice(1):sslice(2)),results)
       ELSE IF (input%bz_integration==BZINT_METHOD_GAUSS) THEN
          CALL fergwt(kpts,input,fmpi,results%neig(:,sslice(1):sslice(2)), results%eig(:,:,sslice(1):sslice(2)),results%ef,results%w_iks(:,:,sslice(1):sslice(2)),results%seigv)
       ELSE IF (input%bz_integration==BZINT_METHOD_TRIA) THEN
          CALL fertri(input,noco,kpts,fmpi%irank, results%neig(:,sslice(1):sslice(2)),nspins,zc,results%eig(:,:,sslice(1):sslice(2)),spindg,&
               results%ef,results%seigv,results%w_iks(:,:,sslice(1):sslice(2)))
       ELSE IF (input%bz_integration==BZINT_METHOD_TETRA) THEN
          CALL fertetra(input,noco,kpts,fmpi,results%neig(:,sslice(1):sslice(2)), results%eig(:,:,sslice(1):sslice(2)),&
                        results%ef,results%w_iks(:,:,sslice(1):sslice(2)),results%seigv)
       ENDIF

       IF (mspin == 2) THEN
          WRITE(oUnit,*) "Different Fermi-energies for both spins:"
          WRITE(oUnit,"(a,f0.3,a,f0.4,a,f0.4,a,f0.4)") "Fixed Moment:" &
               ,input%fixed_moment,"   Difference(EF):",efermi," - ",results%ef,"="&
               ,efermi-results%ef
       ENDIF
       efermi = results%ef
    enddo

    IF (m_spins == 2) nspins = 2

    attributes = ''
    WRITE(attributes(1),'(f20.10)') results%ef
    WRITE(attributes(2),'(a)') 'Htr'
    IF (fmpi%irank.EQ.0) CALL writeXMLElement('FermiEnergy',(/'value','units'/),attributes(1:2))
 ENDIF   

    RETURN
8020 FORMAT (/,'FERMIE:',/,&
         &       10x,'first approx. to ef    (T=0)  :',f11.6,' htr',&
         &       '   (energy of the highest occ. eigenvalue)',/,&
         &       10x,'number of occ. states  (T=0)  :',i11,/,&
         &       10x,'first approx. to seigv (T=0)  :',f11.6,' htr',/,&
         &       10x,'sum of weights of occ. states :',f11.6,/,&
         &       10x,'sum of semicore charge        :',f11.6,' e',/)
  END SUBROUTINE fermie
END MODULE m_fermie
