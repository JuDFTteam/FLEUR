MODULE m_alineso
  USE m_juDFT
  !---------------------------------------------------------------------- 
  ! Set up SO-hamiltonian for 2nd variation (hsohelp and hsoham) and 
  ! diagonalize by lapack routines.
  ! Eigenvalues and vectors (eig_so and zso) are returned 
  !----------------------------------------------------------------------
CONTAINS
  SUBROUTINE alineso(eig_id,&
       mpi,DIMENSION,atoms,sym,&
       input,noco,cell,oneD,&
       rsopp,rsoppd,rsopdp,rsopdpd,nk,&
       rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,&
       usdus,soangl,&
       kveclo,ello,nsize,&
       eig_so,zso)

#include"cpp_double.h"
    USE m_hsohelp
    USE m_hsoham
    USE m_eig66_io, ONLY : read_eig
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_usdus),INTENT(IN)       :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id
    INTEGER, INTENT (IN) :: nk 
    INTEGER, INTENT (OUT):: nsize
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rsopp  (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (IN) :: rsoppd (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (IN) :: rsopdp (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (IN) :: rsopdpd(atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (IN) :: rsoplop (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (IN) :: rsoplopd(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (IN) :: rsopdplo(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (IN) :: rsopplo (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (IN) :: rsoploplop(atoms%ntypd,atoms%nlod,atoms%nlod,2,2)
    COMPLEX, INTENT (IN) :: soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2)
    COMPLEX, INTENT (OUT) :: zso(:,:,:)!(dimension%nbasfcn,2*dimension%neigd,wannierspin)
    REAL,    INTENT (OUT) :: eig_so(2*DIMENSION%neigd),ello(atoms%nlod,atoms%ntypd,DIMENSION%jspd)
    INTEGER, INTENT (OUT) :: kveclo(atoms%nlotot)
    !-odim
    !+odim
    !     ..
    !     .. Local Scalars ..
    TYPE(t_lapw)       :: lapw
    REAL      r2
    INTEGER   i,i1 ,j,jsp,jsp1,k,ne,nn,nn1,nrec,info
    INTEGER   idim_c,idim_r,jsp2,nbas,j1
    CHARACTER vectors 
    LOGICAL   l_file,l_socvec,l_qsgw,l_open
    INTEGER   irec,irecl_qsgw
    COMPLEX   cdum
    !     ..
    !     .. Local Arrays ..
    INTEGER :: nsz(2)
    REAL    :: bkdu(3),eig(DIMENSION%neigd,DIMENSION%jspd),s(3),bkpt(3)
    REAL    :: epar(0:atoms%lmaxd,atoms%ntypd),evac(2)
    REAL,   ALLOCATABLE :: rwork(:)
    COMPLEX,ALLOCATABLE :: cwork(:),chelp(:,:,:,:,:)
    COMPLEX,ALLOCATABLE :: ahelp(:,:,:,:,:),bhelp(:,:,:,:,:)
    COMPLEX,ALLOCATABLE :: zhelp1(:,:),zhelp2(:,:)
    COMPLEX,ALLOCATABLE :: hso(:,:),hsomtx(:,:,:,:)
    COMPLEX,ALLOCATABLE :: sigma_xc_apw(:,:),sigma_xc(:,:)
#ifdef CPP_INVERSION
    REAL,   ALLOCATABLE :: z(:,:,:)
#else
    COMPLEX,ALLOCATABLE :: z(:,:,:)
#endif
    !     ..
    !     .. External Subroutines ..
    EXTERNAL CPP_LAPACK_cheev
    !     ..
    !     .. External Functions ..
    COMPLEX  CPP_BLAS_cdotu,CPP_BLAS_cdotc
    EXTERNAL CPP_BLAS_cdotu,CPP_BLAS_cdotc
    !     ..

    !     read from eigenvalue and -vector file
    !

    INQUIRE (4649,opened=l_socvec)
    INQUIRE (file='fleur.qsgw',exist=l_qsgw)

    ALLOCATE ( z(DIMENSION%nbasfcn,DIMENSION%neigd,DIMENSION%jspd) )

    z(:,:,:)= 0.  
    zso(:,:,:)= CMPLX(0.,0.)

    ALLOCATE(lapw%k1(DIMENSION%nvd,input%jspins))
    ALLOCATE(lapw%k2(DIMENSION%nvd,input%jspins))
    ALLOCATE(lapw%k3(DIMENSION%nvd,input%jspins))
    ALLOCATE(lapw%rk(DIMENSION%nvd,input%jspins))

    DO jsp = 1,input%jspins
       CALL read_eig(&
            eig_id,nk,jsp,&
            bk=bkdu,el=epar,ello=ello(:,:,jsp),&
            evac=evac,neig=ne,eig=eig(:,jsp),&
            nmat=lapw%nmat,nv=lapw%nv(jsp),k1=lapw%k1(:,jsp),k2=lapw%k2(:,jsp),k3=lapw%k3(:,jsp),kveclo=kveclo)
       CALL read_eig(&
            eig_id,nk,jsp,&
            n_start=1,n_end=ne,&
            z=z(:,:ne,jsp))

       ! write(*,*) 'process',irank,' reads ',nk

       bkpt(:) = bkdu(:)
       DO i = 1, lapw%nv(1)
          s(1) = bkpt(1) + lapw%k1(i,1)
          s(2) = bkpt(2) + lapw%k2(i,1)
          s(3) = bkpt(3) + lapw%k3(i,1)
          r2 = DOT_PRODUCT(s,MATMUL(s,cell%bbmat))
          lapw%rk(i,1) = SQRT(r2)
       ENDDO

       IF (ne.GT.DIMENSION%neigd) THEN
          WRITE (6,'(a13,i4,a8,i4)') 'alineso: ne=',ne,' > dimension%neigd=',DIMENSION%neigd
          CALL juDFT_error("alineso: ne > neigd",calledby="alineso")
       ENDIF
       nsz(jsp) = ne
    ENDDO
    !
    ! set up size of e.v. problem in second variation: nsize
    !
    nsize = 0
    DO jsp = 1,input%jspins
       IF (input%jspins.EQ.1) THEN
          nsize = 2*nsz(jsp)
          nsz(2) = nsz(1)
       ELSE
          nsize = nsize + nsz(jsp)
       ENDIF
    ENDDO
    !
    ! set up A and B coefficients
    !
    ALLOCATE ( ahelp(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,atoms%natd,DIMENSION%neigd,DIMENSION%jspd) )
    ALLOCATE ( bhelp(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,atoms%natd,DIMENSION%neigd,DIMENSION%jspd) )
    ALLOCATE ( chelp(-atoms%llod :atoms%llod, DIMENSION%neigd,atoms%nlod,atoms%natd ,DIMENSION%jspd) )
    CALL timestart("alineso SOC: -help") 
    CALL hsohelp(&
         &             DIMENSION,atoms,sym,&
         &             input,lapw,nsz,&
         &             cell,bkpt,&
         &             z,usdus,&
         &             zso,noco,oneD,&
         &             kveclo,&
         &             ahelp,bhelp,chelp)
    CALL timestop("alineso SOC: -help") 
    !
    ! set up hamilton matrix
    !

    CALL timestart("alineso SOC: -ham") 
    ALLOCATE ( hsomtx(2,2,DIMENSION%neigd,DIMENSION%neigd) )
    CALL hsoham(&
         &            atoms,noco,input,nsz,chelp,&
         &            rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,&
         &            ahelp,bhelp,rsopp,rsoppd,rsopdp,rsopdpd,soangl,&
         &            hsomtx)
    DEALLOCATE ( ahelp,bhelp,chelp )
    CALL timestop("alineso SOC: -ham") 
    !
    ! add e.v. on diagonal
    !
    !      write(*,*) '!!!!!!!!!!! remove SOC !!!!!!!!!!!!!!'
    !      hsomtx = 0 !!!!!!!!!!!!
    DO jsp = 1,input%jspins
       DO i = 1,nsz(jsp)
          hsomtx(jsp,jsp,i,i) = hsomtx(jsp,jsp,i,i) +&
               &                           CMPLX(eig(i,jsp),0.)
          IF (input%jspins.EQ.1) THEN
             hsomtx(2,2,i,i) =  hsomtx(2,2,i,i) +&
                  &                           CMPLX(eig(i,jsp),0.)
          ENDIF
       ENDDO
    ENDDO

    !
    !  resort H-matrix 
    !
    ALLOCATE ( hso(2*DIMENSION%neigd,2*DIMENSION%neigd) )
    DO jsp = 1,2
       DO jsp1 = 1,2
          IF (jsp.EQ.1) nn = 0 
          IF (jsp1.EQ.1) nn1 = 0
          IF (jsp.EQ.2) nn = nsz(1)
          IF (jsp1.EQ.2) nn1 = nsz(1)
          !
          DO i = 1,nsz(jsp)
             DO j = 1,nsz(jsp1)
                hso(i+nn,j+nn1) = hsomtx(jsp,jsp1,i,j)
             ENDDO
          ENDDO
          !
       ENDDO
    ENDDO
    DEALLOCATE ( hsomtx )

    !
    !  add Sigma-vxc (QSGW)
    !
    IF( l_qsgw ) THEN
       nbas = lapw%nv(1) + atoms%nlotot
       WRITE(*,'(A,I3,A,I5,A)') 'Read fleur.qsgw  (',nk,',',nbas,')'
       IF( DIMENSION%jspd .EQ. 2 ) STOP 'alineso: GW+noco not implemented.'
       ALLOCATE ( sigma_xc(2*nsz(1),2*nsz(1)) )        
       ALLOCATE ( sigma_xc_apw(nbas,nbas) )
       INQUIRE(667,opened=l_open)
       IF( .NOT.l_open ) THEN
          IF( nk.NE.1 ) STOP 'unit 667 not opened but not at 1st k'
          OPEN(667,file='fleur.qsgw',form='unformatted')
       ELSE IF( nk.EQ.1) THEN
          REWIND(667)
       ENDIF
       DO jsp1 = 1,2
          DO jsp2 = 1,jsp1
             IF(jsp1.EQ.jsp2) THEN
                READ(667) ((sigma_xc_apw(i,j),j=1,i),i=1,nbas)
                DO i = 1,nbas
                   DO j = 1,i-1
                      sigma_xc_apw(j,i) = CONJG(sigma_xc_apw(i,j))
                   ENDDO
                ENDDO
             ELSE
                READ(667) sigma_xc_apw
             ENDIF
             !            write(*,*) 'lo part set to zero!'
             !            sigma_xc_apw(nv+1:,nv+1:) = 0
             i  = nsz(1) * (jsp1-1) + 1 ; i1 = nsz(1) * jsp1
             j  = nsz(1) * (jsp2-1) + 1 ; j1 = nsz(1) * jsp2
#ifdef CPP_INVERSION
             sigma_xc(i:i1,j:j1) = &
                  &        MATMUL (       TRANSPOSE(z(:nbas,:,1))  ,&
                  &        MATMUL ( sigma_xc_apw,   z(:nbas,:,1) ) )
#else
             sigma_xc(i:i1,j:j1) = &
                  &        MATMUL ( CONJG(TRANSPOSE(z(:nbas,:,1))) ,&
                  &        MATMUL ( sigma_xc_apw,   z(:nbas,:,1) ) )
#endif
             hso(i:i1,j:j1) = hso(i:i1,j:j1) + CONJG(sigma_xc(i:i1,j:j1))
             IF(jsp1.NE.jsp2) THEN
                sigma_xc(j:j1,i:i1) = TRANSPOSE(CONJG(sigma_xc(i:i1,j:j1)))
                hso(j:j1,i:i1) = hso(j:j1,i:i1)+CONJG(sigma_xc(j:j1,i:i1))
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE ( sigma_xc_apw )
    ENDIF

    !
    ! diagonalize the hamiltonian using library-routines
    !
    idim_c = 4*DIMENSION%neigd
    idim_r = 6*DIMENSION%neigd

    CALL timestart("alineso SOC: -diag") 

    ALLOCATE ( cwork(idim_c),rwork(idim_r) )

    IF (input%eonly) THEN
       vectors= 'N'
    ELSE
       vectors= 'V'
    ENDIF
    CALL CPP_LAPACK_cheev(vectors,'U',nsize,&
         &                      hso,2*DIMENSION%neigd,&
         &                      eig_so,&
         &                      cwork, idim_c, rwork, &
         &                      info)

    IF (info.NE.0) WRITE (6,FMT=8000) info
8000 FORMAT (' AFTER CPP_LAPACK_cheev: info=',i4)
    CALL timestop("alineso SOC: -diag") 

    DEALLOCATE ( cwork,rwork )

    IF (input%eonly) THEN
       IF(l_socvec)  CALL juDFT_error&
            &        ("EONLY set. Vectors not calculated.",calledby ="alineso")
    ELSE
       ALLOCATE ( zhelp2(DIMENSION%neigd,2*DIMENSION%neigd) )
       !
       ! proj. back to G - space: old eigenvector 'z' to new one 'Z'
       !                                 +
       !  s      ---    s                | z(G,j,s) ...   z(ig,i,jsp)
       ! Z (G) = >     z  (G) * C (i,j)  | Z(G,j,s) ... zso(ig,j,jsp)
       !  j      ---    i                | C(i,j)   ... hso(i ,j)
       !          i                      +
       ! reorder new e.w.  in 2x2 spin space : zhelp(,1),zhelp(,2)
       !
       INQUIRE (file='wann_inp',exist=l_file)

       DO i1 = 1,2

          jsp = i1
          jsp2= i1
          IF (input%jspins.EQ.1) jsp = 1
          IF (input%jspins.EQ.1 .AND..NOT.(l_file.OR.l_socvec)) jsp2=1
          IF (i1.EQ.1) nn = 0
          IF (i1.EQ.2) nn = nsz(1)

          zhelp2(:,:) = 0.d0
          DO j = 1,nsize
             DO i = 1,nsz(jsp)
                zhelp2(i,j) =  CONJG(hso(i+nn,j))
             ENDDO
          ENDDO  ! j

#ifdef CPP_INVERSION
          CALL CPP_BLAS_cgemm("N","N",DIMENSION%nbasfcn,2*dimension%neigd,dimension%neigd,CMPLX(1.d0,0.d0),CMPLX(z(:,:,jsp)),&
               DIMENSION%nbasfcn, zhelp2,DIMENSION%neigd,CMPLX(0.d0,0.d0), zso(1,1,jsp2),DIMENSION%nbasfcn)
#else
          CALL CPP_BLAS_cgemm("N","N",DIMENSION%nbasfcn,2*dimension%neigd,dimension%neigd, CMPLX(1.d0,0.d0),z(:,:,jsp),&
               DIMENSION%nbasfcn, zhelp2,DIMENSION%neigd,CMPLX(0.d0,0.d0), zso(:,:,jsp2),DIMENSION%nbasfcn)
#endif

       ENDDO    !isp

       IF(l_socvec) THEN
          !RS: write SOC vectors to SOCVEC
          WRITE(4649) lapw%nmat,nsize,input%jspins,nsz,2*DIMENSION%neigd,CONJG(hso)
          !CF: write qsgw
          IF(l_qsgw) THEN
             nn = 2*nsz(1)
             sigma_xc = MATMUL ( TRANSPOSE(hso(:nn,:nn)) ,&
                  &                 MATMUL ( sigma_xc , CONJG(hso(:nn,:nn)) ) )
             WRITE(1014) nn
             WRITE(1014) ((sigma_xc(i,j),i=1,j),j=1,nn)
             DEALLOCATE ( sigma_xc )
          ENDIF
       ENDIF

       DEALLOCATE ( zhelp2 )
    ENDIF ! (.NOT.input%eonly)

    DEALLOCATE ( hso,z )
    !
    RETURN

  END SUBROUTINE alineso
END MODULE m_alineso
