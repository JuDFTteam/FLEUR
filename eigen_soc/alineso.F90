MODULE m_alineso
  USE m_juDFT
  !---------------------------------------------------------------------- 
  ! Set up SO-hamiltonian for 2nd variation (hsohelp and hsoham) and 
  ! diagonalize by lapack routines.
  ! Eigenvalues and vectors (eig_so and zso) are returned 
  !----------------------------------------------------------------------
CONTAINS
  SUBROUTINE alineso(eig_id,lapw,&
       mpi,DIMENSION,atoms,sym,kpts,&
       input,noco,cell,oneD, nk, usdus,rsoc,&
       nsize,nmat, eig_so,zso)

#include"cpp_double.h"
    USE m_types
    USE m_hsohelp
    USE m_hsoham
    USE m_eig66_io, ONLY : read_eig
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_lapw),INTENT(IN)        :: lapw
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_usdus),INTENT(IN)       :: usdus
    TYPE(t_rsoc),INTENT(IN)        :: rsoc
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id
    INTEGER, INTENT (IN) :: nk 
    INTEGER, INTENT (OUT):: nsize,nmat
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (OUT) :: zso(:,:,:)!(dimension%nbasfcn,2*dimension%neigd,wannierspin)
    REAL,    INTENT (OUT) :: eig_so(2*DIMENSION%neigd)
    !-odim
    !+odim
    !     ..
    !     .. Local Scalars ..
    REAL      r2
    INTEGER   i,i1 ,j,jsp,jsp1,k,ne,nn,nn1,nrec,info
    INTEGER   idim_c,idim_r,jsp2,nbas,j1,ierr
    CHARACTER vectors 
    LOGICAL   l_socvec,l_qsgw,l_open,l_real
    INTEGER   irec,irecl_qsgw
    INTEGER nat_l, extra, nat_start, nat_stop
    COMPLEX   cdum
    !     ..
    !     .. Local Arrays ..
    INTEGER :: nsz(2)
    REAL    :: eig(DIMENSION%neigd,DIMENSION%jspd),s(3)
    REAL,   ALLOCATABLE :: rwork(:)
    COMPLEX,ALLOCATABLE :: cwork(:),chelp(:,:,:,:,:)
    COMPLEX,ALLOCATABLE :: ahelp(:,:,:,:),bhelp(:,:,:,:)
    COMPLEX,ALLOCATABLE :: zhelp1(:,:),zhelp2(:,:)
    COMPLEX,ALLOCATABLE :: hso(:,:),hsomtx(:,:,:,:)
    COMPLEX,ALLOCATABLE :: sigma_xc_apw(:,:),sigma_xc(:,:)
   
    TYPE(t_mat)::zmat(dimension%jspd)
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

    l_real=sym%invs.and..not.noco%l_noco.and..not.(noco%l_soc.and.atoms%n_u>0)
    zmat%l_real=l_real
    zMat(1:dimension%jspd)%matsize1=lapw%nv(1:dimension%jspd)+atoms%nlotot
    zmat%matsize2=dimension%neigd
   
    INQUIRE (4649,opened=l_socvec)
    INQUIRE (file='fleur.qsgw',exist=l_qsgw)
    if (l_real) THEN
       ALLOCATE (zmat(1)%data_r(zmat(1)%matsize1,DIMENSION%neigd) )
       zmat(1)%data_r(:,:)= 0.  
       if (size(zmat)==2)THEN
          ALLOCATE(zmat(2)%data_r(zmat(2)%matsize1,DIMENSION%neigd) )
          zmat(2)%data_r=0.0
       ENDIF
    else
       ALLOCATE (zmat(1)%data_c(zmat(1)%matsize1,DIMENSION%neigd) )
       zmat(1)%data_c(:,:)= 0.  
       if (size(zmat)==2)THEN
          ALLOCATE(zmat(2)%data_c(zmat(2)%matsize1,DIMENSION%neigd) )
          zmat(2)%data_c=0.0
       ENDIF  
    endif
    zso(:,:,:)= CMPLX(0.,0.)

    DO jsp = 1,input%jspins
       CALL read_eig(&
            eig_id,nk,jsp, neig=ne,eig=eig(:,jsp))
       IF (judft_was_argument("-debugtime")) THEN
          WRITE(6,*) "Non-SOC ev for nk,jsp:",nk,jsp
          WRITE(6,"(6(f10.6,1x))") eig(:ne,jsp)
       ENDIF
       CALL read_eig(&
               eig_id,nk,jsp,&
               n_start=1,n_end=ne,&
               zmat=zmat(jsp))

       ! write(*,*) 'process',irank,' reads ',nk

!!$       DO i = 1, lapw%nv(1)
!!$          s = lapw%bkpt +lapw%gvec(:,i,1)
!!$          r2 = DOT_PRODUCT(s,MATMUL(s,cell%bbmat))
!!$          lapw%rk(i,1) = SQRT(r2)
!!$       ENDDO

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
    ! distribution of (abc)cof over atoms
    !
!
! in case of ev-parallelization, now distribute the atoms:
!
      IF (mpi%n_size > 1) THEN
        nat_l = FLOOR(real(atoms%nat)/mpi%n_size)
        extra = atoms%nat - nat_l*mpi%n_size
        nat_start = mpi%n_rank*nat_l + 1 + extra
        nat_stop  = (mpi%n_rank+1)*nat_l + extra
        IF (mpi%n_rank < extra) THEN
          nat_start = nat_start - (extra - mpi%n_rank)
          nat_stop  = nat_stop - (extra - mpi%n_rank - 1)
        ENDIF
      ELSE
        nat_start = 1
        nat_stop  = atoms%nat
      ENDIF
      nat_l = nat_stop - nat_start + 1
    !
    ! set up A and B coefficients
    !
    ALLOCATE ( ahelp(atoms%lmaxd*(atoms%lmaxd+2),nat_l,DIMENSION%neigd,input%jspins) )
    ALLOCATE ( bhelp(atoms%lmaxd*(atoms%lmaxd+2),nat_l,DIMENSION%neigd,input%jspins) )
    ALLOCATE ( chelp(-atoms%llod :atoms%llod, DIMENSION%neigd,atoms%nlod,nat_l,input%jspins) )
    CALL timestart("alineso SOC: -help") 
    write(*,*) nat_start,nat_stop,nat_l
    CALL hsohelp(&
         &             DIMENSION,atoms,sym,&
         &             input,lapw,nsz,&
         &             cell,&
         &             zmat,usdus,&
         &             zso,noco,oneD,&
         &             nat_start,nat_stop,nat_l,&
         &             ahelp,bhelp,chelp)
    CALL timestop("alineso SOC: -help") 
    !
    ! set up hamilton matrix
    !
    CALL timestart("alineso SOC: -ham") 
#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
    ALLOCATE ( hsomtx(DIMENSION%neigd,DIMENSION%neigd,2,2) )
    CALL hsoham(atoms,noco,input,nsz,dimension%neigd,chelp,rsoc,ahelp,bhelp,&
                nat_start,nat_stop,mpi%n_rank,mpi%n_size,mpi%SUB_COMM,&
                hsomtx)
    write(*,*) 'after hsoham'
    DEALLOCATE ( ahelp,bhelp,chelp )
    CALL timestop("alineso SOC: -ham") 
    IF (mpi%n_rank==0) THEN
    !
    ! add e.v. on diagonal
    !
    !      write(*,*) '!!!!!!!!!!! remove SOC !!!!!!!!!!!!!!'
    !      hsomtx = 0 !!!!!!!!!!!!
    DO jsp = 1,input%jspins
       DO i = 1,nsz(jsp)
          hsomtx(i,i,jsp,jsp) = hsomtx(i,i,jsp,jsp) +&
               &                           CMPLX(eig(i,jsp),0.)
          IF (input%jspins.EQ.1) THEN
             hsomtx(i,i,2,2) =  hsomtx(i,i,2,2) +&
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
          !write(3333,'(2i3,4e15.8)') jsp,jsp1,hsomtx(jsp,jsp1,8,8),hsomtx(jsp,jsp1,32,109) 
          DO i = 1,nsz(jsp)
             DO j = 1,nsz(jsp1)
                hso(i+nn,j+nn1) = hsomtx(i,j,jsp,jsp1)
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
             if (l_real) THEN
                sigma_xc(i:i1,j:j1) = &
                     &        MATMUL (       TRANSPOSE(zmat(1)%data_r(:nbas,:))  ,&
                     &        MATMUL ( sigma_xc_apw,   zmat(1)%data_r(:nbas,:) ) )
else
             sigma_xc(i:i1,j:j1) = &
                  &        MATMUL ( CONJG(TRANSPOSE(zmat(1)%data_c(:nbas,:))) ,&
                  &        MATMUL ( sigma_xc_apw,   zmat(1)%data_c(:nbas,:) ) )
          endif
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

       DO i1 = 1,2

          jsp = i1
          jsp2= i1
          IF (input%jspins.EQ.1) jsp = 1
          IF (input%jspins.EQ.1 .AND..NOT.(input%l_wann.OR.l_socvec)) jsp2=1
          IF (i1.EQ.1) nn = 0
          IF (i1.EQ.2) nn = nsz(1)

          zhelp2(:,:) = 0.d0
          DO j = 1,nsize
             DO i = 1,nsz(jsp)
                zhelp2(i,j) =  CONJG(hso(i+nn,j))
             ENDDO
          ENDDO  ! j

          if (l_real) THEN
             CALL CPP_BLAS_cgemm("N","N",zmat(1)%matsize1,2*dimension%neigd,dimension%neigd,CMPLX(1.d0,0.d0),CMPLX(zmat(jsp)%data_r(:,:)),&
                  zmat(1)%matsize1, zhelp2,DIMENSION%neigd,CMPLX(0.d0,0.d0), zso(1,1,jsp2),zmat(1)%matsize1)
          else
             CALL CPP_BLAS_cgemm("N","N",zmat(1)%matsize1,2*dimension%neigd,dimension%neigd, CMPLX(1.d0,0.d0),zmat(jsp)%data_c(:,:),&
                  zmat(1)%matsize1, zhelp2,DIMENSION%neigd,CMPLX(0.d0,0.d0), zso(:,:,jsp2),zmat(1)%matsize1)
          endif

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

    DEALLOCATE ( hso )
    ENDIF ! (n_rank==0)
    !
    nmat=lapw%nmat
    RETURN

  END SUBROUTINE alineso
END MODULE m_alineso
