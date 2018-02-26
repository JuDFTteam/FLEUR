!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_nonsph
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC hsmt_nonsph
CONTAINS
  SUBROUTINE hsmt_nonsph(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_tlmplm),INTENT(IN)     :: td
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)          :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(IN)            :: chi
    !     .. Array Arguments ..
    REAL,INTENT(IN)               :: fj(:,0:,:),gj(:,0:,:)
    CLASS(t_mat),INTENT(INOUT)     ::hmat
    CALL timestart("non-spherical setup")
    IF (mpi%n_size==1) THEN
       CALL priv_noMPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    ELSE
       CALL priv_MPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    ENDIF
    CALL timestop("non-spherical setup")
  END SUBROUTINE hsmt_nonsph

  SUBROUTINE priv_noMPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(in)   :: chi
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: fj(:,0:,:),gj(:,0:,:)
    CLASS(t_mat),INTENT(INOUT)::hmat

    
    INTEGER:: nn,na,ab_size,l,ll,m
    COMPLEX,ALLOCATABLE:: ab(:,:),ab1(:,:),ab2(:,:)
    real :: rchi

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lmaxd*(atoms%lmaxd+2)+2),ab1(lapw%nv(iintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    IF (iintsp.NE.jintsp) ALLOCATE(ab2(lapw%nv(jintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    IF (hmat%l_real) THEN
       IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
          DEALLOCATE(hmat%data_c)
          ALLOCATE(hmat%data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       ENDIF
       hmat%data_c=0.0
    ENDIF
    
    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
          rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))
          
          CALL hsmt_ab(sym,atoms,noco,isp,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
          !Calculate Hamiltonian
          CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
          !ab1=MATMUL(ab(:lapw%nv(iintsp),:ab_size),td%h_loc(:ab_size,:ab_size,n,isp))
          IF (iintsp==jintsp) THEN
             CALL ZHERK("U","N",lapw%nv(iintsp),ab_size,Rchi,CONJG(ab1),SIZE(ab1,1),1.0,hmat%data_c,SIZE(hmat%data_c,1))
          ELSE  !here the l_ss off-diagonal part starts
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,isp,jintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
             CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab2,SIZE(ab2,1))
             !Multiply for Hamiltonian
             CALL zgemm("N","C",lapw%nv(jintsp),lapw%nv(iintsp),ab_size,chi,ab2,SIZE(ab2,1),ab1,SIZE(ab1,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
          ENDIF
       ENDIF
    END DO
    
    IF (hmat%l_real) THEN
       hmat%data_r=hmat%data_r+REAL(hmat%data_c)
    ENDIF
    
  END SUBROUTINE priv_noMPI


  SUBROUTINE priv_MPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(in)   :: chi
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: fj(:,0:,:),gj(:,0:,:)
    CLASS(t_mat),INTENT(INOUT)::hmat

    
    INTEGER:: nn,na,ab_size,l,ll,m,i,ii
    COMPLEX,ALLOCATABLE:: ab(:,:),ab1(:,:),ab_select1(:,:),ab_select(:,:)
    real :: rchi

<<<<<<< HEAD
                   !
                   !--->             update hamiltonian and overlap matrices
                   nc = 0
                   IF ( noco%l_noco .AND. (n_size>1) ) THEN
                      lapw%nv_tot = lapw%nv(1) + lapw%nv(2)
                   ELSE
                      lapw%nv_tot = lapw%nv(iintsp)
                   ENDIF
                   kii=n_rank
                   DO WHILE(kii<lapw%nv_tot)
                      !DO kii =  n_rank, nv_tot-1, n_size
                      ki = MOD(kii,lapw%nv(iintsp)) + 1
                      bsize=MIN(SIZE(aa_block,1),(lapw%nv(iintsp)-ki)/n_size+1) !Either use maximal blocksize or number of rows left to calculate
                      IF (bsize<1) EXIT !nothing more to do here
                      bsize2=bsize*n_size
                      bsize2=min(bsize2,lapw%nv(iintsp)-ki+1)
                      !aa_block(:bsize,:ki+bsize2-1)=matmul(a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),conjg(transpose(ax(:ki+bsize2-1,0:lmp))))+ &
                      !                              matmul(b(ki:ki+bsize2-1:n_size,0:lmp,iintsp),conjg(transpose(bx(:ki+bsize2-1,0:lmp))))
                      IF (n_size==1) THEN !Make this a special case to avoid copy-in of a array
                         call zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,a(ki,0,iintsp),SIZE(a,1),ax(1,0),SIZE(ax,1),zero,aa_block,SIZE(aa_block,1))
                         call zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,b(ki,0,iintsp),SIZE(a,1),bx(1,0),SIZE(ax,1),one ,aa_block,SIZE(aa_block,1))
                      ELSE
                         CALL zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),SIZE(a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),1),ax(1,0),SIZE(ax,1),zero,aa_block,SIZE(aa_block,1))
                         CALL zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,b(ki:ki+bsize2-1:n_size,0:lmp,iintsp),SIZE(a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),1),bx(1,0),SIZE(ax,1),one,aa_block,SIZE(aa_block,1))
                      ENDIF
                      DO kb=1,bsize
                         IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
                            nc = 1+kii/n_size
                            ii = nc*(nc-1)/2*n_size-(nc-1)*(n_size-n_rank-1)
                            IF ( (n_size==1).OR.(kii+1<=lapw%nv(1)) ) THEN    !
                               aahlp(ii+1:ii+ki) = aahlp(ii+1:ii+ki)+MATMUL(CONJG(ax(:ki,:lmp)),a(ki,:lmp,iintsp))+MATMUL(CONJG(bx(:ki,:lmp)),b(ki,:lmp,iintsp))
                            ELSE                    ! components for <2||2> block unused
                               aa_tmphlp(:ki) = MATMUL(CONJG(ax(:ki,:lmp)),a(ki,:lmp,iintsp))+MATMUL(CONJG(bx(:ki,:lmp)),b(ki,:lmp,iintsp))
                               !--->                   spin-down spin-down part
                               ij = ii + lapw%nv(1)
                               aa_c(ij+1:ij+ki)=aa_c(ij+1:ij+ki)+chi22*aa_tmphlp(:ki)
                               !--->                   spin-down spin-up part, lower triangle
                               ij =  ii
                               aa_c(ij+1:ij+ki)=aa_c(ij+1:ij+ki)+chi21*aa_tmphlp(:ki)
                            ENDIF
                            !-||
                         ELSEIF ( noco%l_noco .AND. noco%l_ss ) THEN
                            IF ( iintsp==1 .AND. jintsp==1 ) THEN
                               !--->                      spin-up spin-up part
                               kjmax = ki
                               chihlp = chi11
                               ii = (ki-1)*(ki)/2
                            ELSEIF ( iintsp==2 .AND. jintsp==2 ) THEN
                               !--->                      spin-down spin-down part
                               kjmax = ki
                               chihlp = chi22
                               ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2+&
                                    lapw%nv(1)+atoms%nlotot
                            ELSE
                               !--->                      spin-down spin-up part
                               kjmax = lapw%nv(1)
                               chihlp = chi21
                               ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2
                            ENDIF
                            aa_c(ii+1:ii+kjmax) = aa_c(ii+1:ii+kjmax) + chihlp*&
                                 (MATMUL(CONJG(ax(:kjmax,:lmp)),a(ki,:lmp,iintsp))+MATMUL(CONJG(bx(:kjmax,:lmp)),b(ki,:lmp,iintsp)))
                         ELSE
                            nc = 1+kii/n_size
                            ii = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                            if (l_real) THEN
                               aa_r(ii+1:ii+ki) = aa_r(ii+1:ii+ki) + aa_block(kb,:ki)
                            ELSE
                               aa_c(ii+1:ii+ki) = aa_c(ii+1:ii+ki) + aa_block(kb,:ki)
                            endif
                            !print*,ii,ki,kb
                            !                           IF (.not.apw(l)) THEN
                            !aa(ii+1:ii+ki) = aa(ii+1:ii+ki) + b(ki,lmp,iintsp)*bx(:ki)
                            !                           ENDIF
                         ENDIF
                         ki=ki+n_size
                         kii=kii+n_size
                      ENDDO
                      !--->             end loop over ki
                   END DO
                   !--->       end loops over interstitial spin
                ENDDO
             ENDDO
          ENDIF              ! atoms%invsat(na) = 0 or 1
          !--->    end loop over equivalent atoms
       END DO
       IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) CALL hsmt_hlptomat(atoms%nlotot,lapw%nv,sub_comm,chi11,chi21,chi22,aahlp,aa_c)
       !---> end loop over atom types (ntype)
    ENDDO ntyploop
=======
    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2),ab1(lapw%nv(iintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2),ab_select(lapw%num_local_cols(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))
>>>>>>> hsmt_simple

    IF (iintsp.NE.jintsp) ALLOCATE(ab_select1(lapw%num_local_cols(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))

    IF (hmat%l_real) THEN
       IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
          DEALLOCATE(hmat%data_c)
          ALLOCATE(hmat%data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       ENDIF
       hmat%data_c=0.0
    ENDIF
    
    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
          rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))
          
          CALL hsmt_ab(sym,atoms,noco,isp,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
          !Calculate Hamiltonian
        
          CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
          IF (iintsp==jintsp) THEN
             !Cut out of ab1 only the needed elements here
             ab_select=ab1(mpi%n_rank+1:lapw%nv(iintsp):mpi%n_size,:)
             CALL zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,CMPLX(1.0,0.0),CONJG(ab1),SIZE(ab1,1),ab_select,lapw%num_local_cols(jintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
          ELSE
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,isp,jintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
             ab_select1=ab(mpi%n_rank+1:lapw%nv(jintsp):mpi%n_size,:)
             CALL zgemm("N","N",lapw%num_local_cols(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab_select1,SIZE(ab_select1,1),td%h_loc(:,:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab_select,SIZE(ab_select,1))
             !Multiply for Hamiltonian
             CALL zgemm("N","t",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,CMPLX(1.0,0.0),CONJG(ab1),SIZE(ab1,1),ab_select,SIZE(ab_select,1)*mpi%n_size,CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))   
          ENDIF
       ENDIF
    END DO
    !delete lower part of matrix
    !i=0
    !DO ii=mpi%n_rank+1,lapw%nv(iintsp),mpi%n_size
    !   i=i+1
    !   hmat%data_c(ii+1:,i)=0.0
    !ENDDO
    IF (hmat%l_real) THEN
       hmat%data_r=hmat%data_r+REAL(hmat%data_c)
    ENDIF
    
  END SUBROUTINE priv_MPI

  
END MODULE m_hsmt_nonsph
