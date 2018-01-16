!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_lapw
  USE m_types_misc
  USE m_judft
  PRIVATE
  TYPE t_lapw
     INTEGER :: nv(2),num_local_cols(2)
     INTEGER :: nv_tot
     INTEGER :: nmat
     INTEGER :: nlotot
     INTEGER,ALLOCATABLE:: k1(:,:)
     INTEGER,ALLOCATABLE:: k2(:,:)
     INTEGER,ALLOCATABLE:: k3(:,:)
     INTEGER,ALLOCATABLE:: gvec(:,:,:) !replaces k1,k2,k3
     INTEGER,ALLOCATABLE:: kp(:,:)
     REAL,ALLOCATABLE::rk(:,:)
     REAL,ALLOCATABLE::gk(:,:,:)
     REAL,ALLOCATABLE::vk(:,:,:)
     INTEGER,ALLOCATABLE::matind(:,:)
     INTEGER,ALLOCATABLE::index_lo(:,:)
     INTEGER,ALLOCATABLE::kvec(:,:,:)
     INTEGER,ALLOCATABLE::kveclo(:)
     REAL   :: bkpt(3)
   CONTAINS
     PROCEDURE,PASS :: init =>lapw_init
     PROCEDURE,PASS :: alloc =>lapw_alloc
     PROCEDURE,PASS :: phase_factors =>lapw_phase_factors
  END TYPE t_lapw
  PUBLIC :: t_lapw


CONTAINS
  SUBROUTINE lapw_alloc(lapw,cell,input,noco)
    !
    !*********************************************************************
    !     determines dimensions of the lapw basis set with |k+G|<rkmax.
    !     bkpt is the k-point given in internal units
    !*********************************************************************
    USE m_boxdim
    IMPLICIT NONE
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_noco),INTENT(IN)      :: noco
    CLASS(t_lapw),INTENT(INOUT)  :: lapw


    INTEGER j1,j2,j3,mk1,mk2,mk3,nv
    INTEGER ispin,nvh(2)

    REAL arltv1,arltv2,arltv3,rkm,rk2,r2,s(3)
    ! ..
    !
    !------->          ABBREVIATIONS
    !

    !   rkmax       : cut-off for |g+k|
    !   arltv(i)    : length of reciprical lattice vector along
    !                 direction (i)
    !
    !---> Determine rkmax box of size mk1, mk2, mk3,
    !     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
    !
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

    !     (add 1+1 due to integer rounding, strange k_vector in BZ)
    mk1 = int(input%rkmax/arltv1) + 2
    mk2 = int(input%rkmax/arltv2) + 2
    mk3 = int(input%rkmax/arltv3) + 2

    rkm = input%rkmax
    rk2 = rkm*rkm
    !---> obtain vectors
    !---> in a spin-spiral calculation different basis sets are used for
    !---> the two spin directions, because the cutoff radius is defined
    !---> by |G + k +/- qss/2| < rkmax.
    nvh(2)=0
    DO ispin = 1,MERGE(2,1,noco%l_ss)
       nv = 0
       DO j1 = -mk1,mk1
          DO j2 = -mk2,mk2
             DO j3 = -mk3,mk3
                s = lapw%bkpt + (/j1,j2,j3/) + (2*ispin - 3)/2.0*noco%qss
                r2 = dot_PRODUCT(MATMUL(s,cell%bbmat),s)
                IF (r2.LE.rk2)  nv = nv + 1
             END DO
          END DO
       END DO
       nvh(ispin)  = nv
    END DO
    nv  = MAX(nvh(1),nvh(2))

    IF (ALLOCATED(lapw%rk)) THEN
       IF (SIZE(lapw%rk)==nv) THEN
          RETURN !
       ELSE
          DEALLOCATE(lapw%rk,lapw%gvec,lapw%vk,lapw%gk,lapw%matind)
       ENDIF
    ENDIF
    ALLOCATE(lapw%rk(nv,input%jspins) )
    ALLOCATE(lapw%gvec(3,nv,input%jspins))
    ALLOCATE(lapw%vk(3,nv,input%jspins))
    ALLOCATE(lapw%gk(3,nv,input%jspins))
    ALLOCATE(lapw%k1(nv,2)) !shpuld be removed
    ALLOCATE(lapw%k2(nv,2)) !
    ALLOCATE(lapw%k3(nv,2)) !
    ALLOCATE(lapw%matind(nv,2))

    lapw%rk = 0 ; lapw%gvec = 0 ;lapw%nv=0
  END SUBROUTINE lapw_alloc


  SUBROUTINE lapw_init(lapw,input,noco,kpts,atoms,sym,&
       nk,cell,l_zref,mpi)
    USE m_types_mpi
    USE m_sort
    USE m_boxdim
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_mpi),INTENT(IN),OPTIONAL:: mpi
    CLASS(t_lapw),INTENT(INOUT)    :: lapw
    !     .. 
    !     .. Scalar Arguments ..
    INTEGER, INTENT  (IN) :: nk
    LOGICAL, INTENT (IN)  :: l_zref
    !     ..
    !     .. Array Arguments ..
    !     ..
    !     .. Local Scalars ..
    REAL arltv1,arltv2,arltv3,r2,rk2,rkm,r2q,gla,eps
    INTEGER i,j,j1,j2,j3,k,l ,mk1,mk2,mk3,n,ispin,gmi,m,nred,n_inner,n_bound
    !     ..
    !     .. Local Arrays ..
    REAL                :: s(3),sq(3)
    REAL,ALLOCATABLE    :: rk(:),rkq(:)
    INTEGER,ALLOCATABLE :: gvec(:,:),index3(:)


    !     ..
    !---> in a spin-spiral calculation different basis sets are used for
    !---> the two spin directions, because the cutoff radius is defined
    !---> by |G + k +/- qss/2| < rkmax.

    IF (nk>kpts%nkpt) THEN
       lapw%bkpt(:)=kpts%bkf(:,nk)
    ELSE
       lapw%bkpt(:) = kpts%bk(:,nk)
    ENDIF

    CALL lapw%alloc(cell,input,noco)

    ALLOCATE(gvec(3,SIZE(lapw%gvec,2)))
    ALLOCATE(rk(SIZE(lapw%gvec,2)),rkq(SIZE(lapw%gvec,2)))
    ALLOCATE(index3(SIZE(lapw%gvec,2)))


    !---> Determine rkmax box of size mk1, mk2, mk3,
    !     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
    !     arltv(i) length of reciprical lattice vector along direction (i)
    !
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

    !     (add 1+1 due to integer rounding, strange k_vector in BZ)
    mk1 = int( input%rkmax/arltv1 ) + 4
    mk2 = int( input%rkmax/arltv2 ) + 4
    mk3 = int( input%rkmax/arltv3 ) + 4

    rk2 = input%rkmax*input%rkmax
    !---> if too many basis functions, reduce rkmax
    spinloop:DO ispin = 1,input%jspins
       rk2 = rkm*rkm
       !--->    obtain vectors
       n = 0
       DO  j1 = -mk1,mk1
          DO  j2 = -mk2,mk2
             DO  j3 = -mk3,mk3
                s=lapw%bkpt+(/j1,j2,j3/)+(2*ispin - 3)/2.0*noco%qss
                sq = lapw%bkpt+ (/j1,j2,j3/)
                r2 = dot_PRODUCT(s,MATMUL(s,cell%bbmat))
                r2q = dot_PRODUCT(sq,MATMUL(sq,cell%bbmat))
                IF (r2.LE.rk2) THEN
                   n = n + 1
                   gvec(:,n) = (/j1,j2,j3/)
                   rk(n) = SQRT(r2)
                   rkq(n) = SQRT(r2q)
                END IF
             ENDDO
          ENDDO
       ENDDO
       lapw%nv(ispin) = n

       !Sort according to k+g
       CALL sort(lapw%nv(ispin),rkq,index3)
       DO n=1,lapw%nv(ispin)
          lapw%gvec(:,n,ispin) = gvec(:,index3(n))
          lapw%rk(n,ispin)     = rk(index3(n))
       ENDDO

       !+gu
       !--->    determine pairs of K-vectors, where K_z = K'_-z to use 
       !--->    z-reflection
       IF (l_zref) THEN
          n=0
          DO i=1,lapw%nv(ispin)
             DO j=1,i
                IF (ALL(lapw%gvec(1:2,i,ispin).EQ.lapw%gvec(1:2,j,ispin)).AND.&
                     (lapw%gvec(3,i,ispin).EQ.-lapw%gvec(3,j,ispin))) THEN
                   n=n+1 
                   lapw%matind(n,1)=i
                   lapw%matind(n,2)=j
                ENDIF
             ENDDO
          ENDDO
          nred=n
          IF (PRESENT(mpi)) THEN
             IF (mpi%n_size.GT.1) THEN
             !
             !--->     order K's in sequence K_1,...K_n | K_0,... | K_-1,....K_-n
             !
             n_inner = lapw%nv(ispin) - nred
             IF (MOD(nred,mpi%n_size).EQ.0) THEN
                n_bound = nred
             ELSE
                n_bound = (1+INT( nred/mpi%n_size ))*mpi%n_size
             ENDIF
             IF (lapw%nv(ispin) - nred + n_bound.GT.SIZE(lapw%gvec,2)) THEN
                CALL juDFT_error("BUG:z-ref & ev || : dimension too small!" ,calledby ="types_lapw")
             ENDIF

             i = 1
             j = 1
             DO n = 1, nred 
                IF (lapw%matind(n,1).EQ.lapw%matind(n,2)) THEN
                   index3(lapw%matind(n,1)) = n_inner + i
                   i = i + 1
                ELSE
                   index3(lapw%matind(n,1)) = j
                   index3(lapw%matind(n,2)) = j + n_bound
                   j = j + 1
                ENDIF
             ENDDO
             !--->          resort the rk,k1,k2,k3 and lapw%matind arrays:
             DO n = 1, lapw%nv(ispin)
                rk(n)  = lapw%rk(n,ispin)
                gvec(:,n) = lapw%gvec(:,n,ispin)
             ENDDO
             DO n = lapw%nv(ispin), 1, -1
                lapw%rk(index3(n),ispin) = rk(n)
                lapw%gvec(:,index3(n),ispin) = gvec(:,n)
             ENDDO
             DO n = nred + 1, n_bound
                lapw%rk(n,ispin) = lapw%rk(lapw%nv(ispin),ispin)
                lapw%gvec(:,n,ispin) = lapw%gvec(:,lapw%nv(ispin),ispin)
             ENDDO
             lapw%nv(ispin) = lapw%nv(ispin) - nred + n_bound
          ENDIF
       ENDIF
       ENDIF

       IF (noco%l_ss) THEN  ! sort additionally like in strgn1... gb
          i = 1
          gla = 0.
          rk(1) = 0.0
          eps=1.e-10
          DO  k = 1,lapw%nv(ispin)
             IF (rkq(k)-gla.GE.eps) i=i+1
             gla = rkq(k)
             gmi = (mk1+lapw%gvec(1,k,ispin)) + (mk2+lapw%gvec(2,k,ispin))*(2*mk1+1) +&
                  (mk3+lapw%gvec(3,k,ispin))*(2*mk1+1)*(2*mk2+1)
             rk(k) = i * (9.+(2*mk1+1)*(2*mk2+1)*(2*mk3+1)) + gmi
          ENDDO
          CALL sort(lapw%nv(ispin),rk,index3)
          DO  k = 1,lapw%nv(ispin)
             gvec(:,k) = lapw%gvec(:,index3(k),ispin)
             rk(k) =  lapw%rk(index3(k),ispin)
          ENDDO
          lapw%gvec(:,:lapw%nv(ispin),ispin) = gvec(:,:lapw%nv(ispin))
          lapw%rk(:lapw%nv(ispin),ispin) = rk(:lapw%nv(ispin))
       ENDIF
       !-gu
       DO k=1,lapw%nv(ispin)
          lapw%vk(:,k,ispin)=lapw%bkpt+lapw%gvec(:,k,ispin)+(ispin-1.5)*noco%qss
          lapw%gk(:,k,ispin)=MATMUL(TRANSPOSE(cell%bmat),lapw%vk(:,k,ispin))/MAX (lapw%rk(k,ispin),1.0e-30)
       ENDDO

       IF (.NOT.noco%l_ss.AND.input%jspins==2) THEN
          !Second spin is the same
          lapw%nv(2)=lapw%nv(1)
          lapw%gvec(:,:,2)=lapw%gvec(:,:,1)
          lapw%rk(:,2)=lapw%rk(:,1)
          lapw%vk(:,:,2)=lapw%vk(:,:,1)
          lapw%gk(:,:,2)=lapw%gk(:,:,1)
          EXIT spinloop
       END IF

    ENDDO spinloop
    !should be removed later...
    lapw%k1=lapw%gvec(1,:,:)
    lapw%k2=lapw%gvec(2,:,:)
    lapw%k3=lapw%gvec(3,:,:)

    lapw%num_local_cols=lapw%nv+atoms%nlotot
    !TODO needs to be adjusted in MPI case

    IF (ANY(atoms%nlo>0)) CALL priv_lo_basis_setup(lapw,atoms,sym,noco,cell)

  CONTAINS

  SUBROUTINE priv_lo_basis_setup(lapw,atoms,sym,noco,cell)
    USE m_vecforlo
    IMPLICIT NONE
    TYPE(t_lapw),INTENT(INOUT):: lapw
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_noco),INTENT(IN)   :: noco


    INTEGER:: n,na,nn,np,lo,nkvec_sv,nkvec(atoms%nlod,2),iindex
    IF (.NOT.ALLOCATED(lapw%kvec)) THEN
       ALLOCATE(lapw%kvec(2*(2*atoms%llod+1),atoms%nlod,atoms%ntype))
       ALLOCATE(lapw%index_lo(atoms%nlod,atoms%ntype))
       ALLOCATE(lapw%kveclo(atoms%nlotot))
    ENDIF
    iindex=0
    na=0
    nkvec_sv=0
    DO n=1,atoms%ntype
       DO nn=1,atoms%neq(n)
          na=na+1
          !np = MERGE(oneD%ods%ngopr(na),sym%invtab(atoms%ngopr(na)),oneD%odi%d1)
          np=sym%invtab(atoms%ngopr(na))
          CALL priv_vec_for_lo(atoms,sym,na,n,np,noco,lapw,cell,nkvec)
          DO lo = 1,atoms%nlo(n)
             lapw%index_lo(lo,n)=iindex
             iindex=iindex+nkvec(lo,1)
             lapw%kveclo(iindex+1:iindex+nkvec(lo,1)) = lapw%kvec(1:nkvec(lo,1),lo,n)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE priv_lo_basis_setup
    
  END SUBROUTINE lapw_init


  SUBROUTINE lapw_phase_factors(lapw,iintsp,tau,qss,cph)
    USE m_constants
    IMPLICIT NONE
    CLASS(t_lapw),INTENT(in):: lapw
    INTEGER,INTENT(IN)     :: iintsp
    REAL,INTENT(in)        :: tau(3),qss(3)
    COMPLEX,INTENT(out)    :: cph(:)

    INTEGER:: k
    REAL:: th
    DO k = 1,lapw%nv(iintsp)
       th= DOT_PRODUCT(lapw%gvec(:,k,iintsp)+(iintsp-1.5)*qss,tau)
       cph(k) = CMPLX(COS(tpi_const*th),-SIN(tpi_const*th))
    END DO
  END SUBROUTINE lapw_phase_factors

  
  SUBROUTINE priv_vec_for_lo(atoms,sym,na,&
       n,np,noco, lapw,cell, nkvec)
    USE m_constants,ONLY: tpi_const,fpi_const
    USE m_orthoglo
    USE m_ylm
    USE m_types_misc
    IMPLICIT NONE
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(INOUT):: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,n,np 
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (INOUT):: nkvec(atoms%nlod,2)
    !     ..
    !     .. Local Scalars ..
    COMPLEX term1 
    REAL th,con1
    INTEGER l,lo ,mind,ll1,lm,iintsp,k,nkmin,ntyp,lmp,m,nintsp
    LOGICAL linind,enough,l_lo1
    !     ..
    !     .. Local Arrays ..
    REAL qssbti(3),bmrot(3,3),v(3),vmult(3)
    REAL :: gkrot(3,SIZE(lapw%gk,2),2)
    REAL :: rph(SIZE(lapw%gk,2),2)
    REAL :: cph(SIZE(lapw%gk,2),2)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    COMPLEX cwork(-2*atoms%llod:2*atoms%llod+1,2*(2*atoms%llod+1),atoms%nlod ,2)
    !     ..
    !     .. Data statements ..
    REAL, PARAMETER :: eps = 1.0E-30
    REAL, PARAMETER :: linindq = 1.0e-4

    con1=fpi_const/SQRT(cell%omtil)
    ntyp = n
    nintsp=MERGE(2,1,noco%l_ss)
    DO iintsp = 1,nintsp
       IF (iintsp.EQ.1) THEN
          qssbti = - noco%qss/2
       ELSE
          qssbti = + noco%qss/2
       ENDIF

       !--->    set up phase factors
       DO k = 1,lapw%nv(iintsp)
          th= tpi_const*DOT_PRODUCT((/lapw%k1(k,iintsp),lapw%k2(k,iintsp),lapw%k3(k,iintsp)/)+qssbti,atoms%taual(:,na))
          rph(k,iintsp) = COS(th)
          cph(k,iintsp) = -SIN(th)
       END DO

       IF (np.EQ.1) THEN
          gkrot(:,:,iintsp)=lapw%gk(:,:,iintsp)
       ELSE
          bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
          DO k = 1,lapw%nv(iintsp)
             !-->           apply the rotation that brings this atom into the
             !-->           representative (this is the definition of ngopr(na))
             !-->           and transform to cartesian coordinates
             v(:) = lapw%vk(:,k,iintsp)
             gkrot(:,k,iintsp) = MATMUL(v,bmrot)
          END DO
       END IF
       !--->   end loop over interstitial spin
    ENDDO

    nkvec(:,:) = 0
    cwork(:,:,:,:) = CMPLX(0.0,0.0)
    enough=.FALSE.
    DO k = 1,MIN(lapw%nv(1),lapw%nv(nintsp))
       IF (ANY(lapw%rk(k,:nintsp).LT.eps)) CYCLE
       IF (.NOT.enough) THEN
          DO iintsp = 1,nintsp

             !-->        generate spherical harmonics
             vmult(:) =  gkrot(:,k,iintsp)
             CALL ylm4(atoms%lnonsph(ntyp),vmult, ylm)
                enough = .TRUE.
                term1 = con1* ((atoms%rmt(ntyp)**2)/2)* CMPLX(rph(k,iintsp),cph(k,iintsp))
                DO lo = 1,atoms%nlo(ntyp)
                   IF (atoms%invsat(na).EQ.0) THEN
                      IF ((nkvec(lo,iintsp)).LT. (2*atoms%llo(lo,ntyp)+1)) THEN
                         enough = .FALSE.
                         nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                         l = atoms%llo(lo,ntyp)
                         ll1 = l*(l+1) + 1
                         DO m = -l,l
                            lm = ll1 + m
                            cwork(m,nkvec(lo,iintsp),lo,iintsp) = term1*ylm(lm)
                         END DO
                         CALL orthoglo(&
                              sym%invs,atoms,nkvec(lo,iintsp),lo,l,linindq,.FALSE., cwork(-2*atoms%llod,1,1,iintsp),linind)
                         IF (linind) THEN
                            lapw%kvec(nkvec(lo,iintsp),lo,ntyp) = k
                         ELSE
                            nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                         ENDIF
                      ENDIF
                   ELSE
                      IF ((atoms%invsat(na).EQ.1) .OR. (atoms%invsat(na).EQ.2)) THEN
                         IF (nkvec(lo,iintsp).LT.2*(2*atoms%llo(lo,ntyp)+1)) THEN
                            enough = .FALSE.
                            nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                            l = atoms%llo(lo,ntyp)
                            ll1 = l*(l+1) + 1
                            DO m = -l,l
                               lm = ll1 + m
                               mind = -l + m
                               cwork(mind,nkvec(lo,iintsp),lo,iintsp) = term1*ylm(lm)
                               mind = l + 1 + m
                               lmp = ll1 - m
                               cwork(mind,nkvec(lo,iintsp),lo,iintsp) = ((-1)** (l+m))*CONJG(term1*ylm(lmp))
                            END DO
                            CALL orthoglo(&
                                 sym%invs,atoms,nkvec(lo,iintsp),lo,l,linindq,.TRUE., cwork(-2*atoms%llod,1,1,iintsp),linind)
                            IF (linind) THEN
                               lapw%kvec(nkvec(lo,iintsp),lo,ntyp) = k
                               !                          write(*,*) nkvec(lo,iintsp),k,' <- '
                            ELSE
                               nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                            END IF
                         END IF
                      END IF
                   END IF
                END DO
                IF ((k.EQ.lapw%nv(iintsp)) .AND. (.NOT.enough)) THEN
                   WRITE (6,FMT=*) 'vec_for_lo did not find enough linearly independent'
                   WRITE (6,FMT=*) 'clo coefficient-vectors. the linear independence'
                   WRITE (6,FMT=*) 'quality, linindq, is set: ',linindq
                   WRITE (6,FMT=*) 'this value might be to large.'
                   WRITE(*,*) na,k,lapw%nv 
                   CALL juDFT_error("not enough lin. indep. clo-vectors" ,calledby ="vec_for_lo")
                END IF
             ! -- >        end of abccoflo-part           
          ENDDO
       ENDIF

       ! -->    check whether we have already enough k-vecs
       enough=.TRUE.
       DO lo = 1,atoms%nlo(ntyp)
          IF (nkvec(lo,1).EQ.nkvec(lo,nintsp)) THEN   ! k-vec accepted by both spin channels
             IF (atoms%invsat(na).EQ.0) THEN
                IF ( nkvec(lo,1).LT.(2*atoms%llo(lo,ntyp)+1) ) THEN 
                   enough=.FALSE.
                ENDIF
             ELSE
                IF ( nkvec(lo,1).LT.(2*(2*atoms%llo(lo,ntyp)+1))  ) THEN
                   enough=.FALSE.
                ENDIF
             ENDIF
          ELSE
             nkmin = MIN(nkvec(lo,1),nkvec(lo,nintsp)) ! try another k-vec
             nkvec(lo,1) = nkmin ; nkvec(lo,nintsp) = nkmin
             enough=.FALSE.
          ENDIF
       ENDDO
       IF ( enough ) THEN
          !           DO  lo = 1,nlo(ntyp)
          !             DO l = 1, nkvec(lo,1)
          !              write(*,*) lo,l,kvec(l,lo)
          !             ENDDO
          !           ENDDO
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE priv_vec_for_lo
END MODULE m_types_lapw
